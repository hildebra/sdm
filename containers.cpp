/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand
email: Falk.Hildebrand@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <algorithm>
#include <cassert>
#include <deque>
#include <future>
#include <limits>
#include <filesystem>
#include <sstream>
#include "OutputFileStream.h"
#include <mutex>
#include "Common.h"
#include "containers.h"
#include "PartitionWorkPool.h"
#include "ThreadPool.h"
#include "CoarseDereplication.h"
#include "Filters.h"
#include "OutputStreamer.h"

using namespace std;

namespace {
	OptContainer* g_opt_container = nullptr;

	// Runs task(0 .. count-1) on the pool and waits for all of them; the first
	// failure is rethrown once every submitted task has finished. Never call this
	// from one of the pool's own workers.
	template <class Task>
	void runOnPool(ThreadPool& pool, size_t count, const Task& task) {
		std::vector<std::future<void>> done;
		done.reserve(count);
		std::exception_ptr failure;
		try {
			for (size_t i = 0; i < count; ++i) { done.push_back(pool.submit([&task, i]() { task(i); })); }
		} catch (...) { failure = std::current_exception(); }
		for (auto& future : done) {
			try { future.get(); } catch (...) { if (!failure) { failure = std::current_exception(); } }
		}
		if (failure) { std::rethrow_exception(failure); }
	}

	// Sorts `items` in place by `less`. With a pool and enough items, every
	// segment is split at its median (std::nth_element), all segments of a level
	// in parallel, until there is one segment per worker; the segments are then
	// sorted in parallel. Needs no memory per item. `less` must be a strict total
	// order (no equivalent items), so the result is exactly std::sort's.
	template <class T, class Less>
	void sortOnPool(std::vector<T>& items, const Less& less, ThreadPool* pool) {
		constexpr size_t kSegment = 16384;
		if (!pool || items.size() < 2 * kSegment) { std::sort(items.begin(), items.end(), less); return; }
		std::vector<size_t> bounds{0, items.size()};
		while (bounds.size() - 1 < pool->size() && items.size() / (bounds.size() - 1) >= 2 * kSegment) {
			const size_t segments = bounds.size() - 1;
			std::vector<size_t> next(2 * segments + 1);
			for (size_t s = 0; s < segments; ++s) {
				next[2 * s] = bounds[s];
				next[2 * s + 1] = bounds[s] + (bounds[s + 1] - bounds[s]) / 2;
			}
			next[2 * segments] = items.size();
			runOnPool(*pool, segments, [&](size_t s) {
				std::nth_element(items.begin() + next[2 * s], items.begin() + next[2 * s + 1], items.begin() + next[2 * s + 2], less);
			});
			bounds.swap(next);
		}
		runOnPool(*pool, bounds.size() - 1, [&](size_t s) { std::sort(items.begin() + bounds[s], items.begin() + bounds[s + 1], less); });
	}

	// Formatted output of one dereplication output chunk, kept in fixed-size
	// segments taken from a shared pool. Growing and freeing a string per
	// chunk on many threads leaves large freed buffers in each worker's malloc
	// arena; recycled segments keep the memory at what is actually in flight.
	class SegmentPool {
	public:
		static constexpr size_t kSegmentSize = size_t(1) << 20;
		using Segment = std::unique_ptr<char[]>;
		Segment take() {
			std::lock_guard<std::mutex> lock(mutex_);
			if (free_.empty()) { return Segment(new char[kSegmentSize]); }
			Segment segment = std::move(free_.back());
			free_.pop_back();
			return segment;
		}
		void give(Segment segment) {
			std::lock_guard<std::mutex> lock(mutex_);
			free_.push_back(std::move(segment));
		}
	private:
		std::mutex mutex_;
		std::vector<Segment> free_;
	};

	class SegmentBuffer final : public std::streambuf {
	public:
		explicit SegmentBuffer(SegmentPool& pool) : pool_(pool) {}
		SegmentBuffer(const SegmentBuffer&) = delete;
		SegmentBuffer& operator=(const SegmentBuffer&) = delete;
		~SegmentBuffer() override { release(); }
		// Writes the buffered bytes and returns the segments to the pool.
		void writeTo(std::ostream& out) {
			for (size_t i = 0; i < segments_.size(); ++i) {
				const size_t used = i + 1 == segments_.size() ? static_cast<size_t>(pptr() - pbase()) : SegmentPool::kSegmentSize;
				out.write(segments_[i].get(), static_cast<std::streamsize>(used));
			}
			release();
		}
	protected:
		int_type overflow(int_type ch) override {
			if (traits_type::eq_int_type(ch, traits_type::eof())) { return traits_type::not_eof(ch); }
			grow();
			*pptr() = traits_type::to_char_type(ch);
			pbump(1);
			return ch;
		}
		std::streamsize xsputn(const char* data, std::streamsize count) override {
			std::streamsize written = 0;
			while (written < count) {
				if (pptr() == epptr()) { grow(); }
				const std::streamsize amount = std::min<std::streamsize>(epptr() - pptr(), count - written);
				std::memcpy(pptr(), data + written, static_cast<size_t>(amount));
				pbump(static_cast<int>(amount));
				written += amount;
			}
			return written;
		}
	private:
		void grow() {
			segments_.push_back(pool_.take());
			char* begin = segments_.back().get();
			setp(begin, begin + SegmentPool::kSegmentSize);
		}
		void release() {
			for (auto& segment : segments_) { pool_.give(std::move(segment)); }
			segments_.clear();
			setp(nullptr, nullptr);
		}
		SegmentPool& pool_;
		std::vector<SegmentPool::Segment> segments_;
	};

	struct SegmentStream {
		explicit SegmentStream(SegmentPool& pool) : buffer(pool), stream(&buffer) {}
		SegmentBuffer buffer;
		std::ostream stream;
	};
}

void setGlobalOptContainer(OptContainer* options) {
	g_opt_container = options;
}

OptContainer* getGlobalOptContainer() {
	return g_opt_container;
}

OptContainer& globalOptContainer() {
	if (g_opt_container == nullptr) {
		throw std::runtime_error("Global option container is not initialized.");
	}
	return *g_opt_container;
}





// trim and is_digits are provided by Common.cpp


//



ReadSubset::ReadSubset(const string inf, const string default_outfile):
RemainderStrPos(-1), newHD(0), outFiles(0), outFilesIdx(0) {
	string line;
	ifstream in(inf.c_str());
	if (!in){
		cerr << "Could not find " << inf << " read subset file. Exiting.\n"; exit(90);
	}
    int ini_ColPerRow(0), cnt(0), skips(0);

	//check read subset format
	while (!safeGetline(in, line).eof()) {
		// Same rows as the parser below reads: it skips comments and lines under 5 characters.
		if (line.substr(0, 1) == "#" || line.length() < 5){ skips++; continue; }
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss, segments, '\t')) {
			ColsPerRow++;
		}

		if (cnt == 0){
			ini_ColPerRow = ColsPerRow;
		}
		else {
			if (ColsPerRow != ini_ColPerRow){
				cerr << "Number of columns in read subset file on line " << cnt + skips + 1 << " is " << ColsPerRow << ". Expected " << ini_ColPerRow << " columns.\n";
				exit(91);
			}
			// Every row is checked: a short row after the first 1000 misaligned
			// the header and output-file columns (out-of-bounds outFilesIdx).
		}
		cnt++;
	}
	if (ini_ColPerRow == 0){
		cerr << "Read Subset File exists, but appears to be badly formated (0 columns detected). Exiting\n"; exit(92);
	}
	if (cnt == 0){
		cerr << "Read Subset File exists, but appears to be badly formated (0 lines detected). Exiting\n"; exit(92);
	}
	in.clear();
	in.seekg(0, ios::beg);
	//extract read subset content
	cnt = 0;
	map<string, int> tmpFiles;
	map<string, int>::iterator tmpFilIT;
	unordered_map <string, int>::iterator TarIT;
	//parameters were set for in matrix, now read line by line
	while (!safeGetline(in, line).eof()) {
		//	while(getline(in,line,'\n')) {
		if (line.substr(0, 1) == "#"){ continue; } // comment/header, as in the format check above
		if (line.length() < 5){ continue; }
		stringstream ss; string segments; int cnt2 = 0;
		ss << line;
		while (getline(ss, segments, '\t')) {
			if (cnt2 == 0){
				//free target RD from paired end info
				remove_paired_info(segments);
				TarIT = Targets.find(segments);
				if (TarIT == Targets.end()){
					Targets[segments] = cnt;
				} else {
					break;
				}
			}
			else if (cnt2 == 1){
				newHD.push_back(segments);
			}
			else if (cnt2 == 2){
				uint Idx(0);
				//1: test if outfile already exists
				tmpFilIT = tmpFiles.find(segments);
				if (tmpFilIT != tmpFiles.end()){
					Idx = (*tmpFilIT).second;
				}
				else {//create new key
					Idx = (uint) outFiles.size();
					outFiles.push_back(segments);
					tmpFiles[segments] = Idx;
				}
				//2: add the index to outfile to this position
				outFilesIdx.push_back(Idx);
			}
			cnt2++;
		}
		if (cnt2 == 1) {//no specific header
			newHD.push_back("");
			tmpFilIT = tmpFiles.find("Default");
			uint Idx(0);
			if (tmpFilIT != tmpFiles.end()) {
				Idx = (*tmpFilIT).second;
			} else {//create new key
				Idx = (uint)outFiles.size();
				outFiles.push_back("Default");
				tmpFiles["Default"] = Idx;
			}

			outFilesIdx.push_back(Idx);
		}
		if (cnt2>0) { cnt++; }
	}
	//is extra file info in read subset file?
	if (ini_ColPerRow < 3){
		outFilesIdx.resize(Targets.size(), 0);
		outFiles.resize(1, default_outfile);
	}

}

void ReadSubset::findMatches(shared_ptr<InputStreamer> IS, OutputStreamer* MD, bool mocatFix) {
	vector<shared_ptr<DNA>> match; //shared_ptr<DNA> match2(NULL);
	bool cont(true), cont2(true);
	int idx(0);
	int pairs = IS->pairNum();
	unordered_map <string, int>::iterator SEEK;
	bool b_doHD = newHD.size() > 0;
	//bool sync(false);//meaningless placeholder
	while (cont) {
		match = IS->getDNApairs();
		if (match[0] == NULL) { cont = false; break; }
		//if (pairs > 1) {			match2 = IS->getDNA(cont2, 1);		}
		string curID = match[0]->getPositionFreeId();
		if (mocatFix && curID.length()>5 ) {
			curID = header_stem(curID);
			//cerr << curID << endl;
			//exit(0);
		} 
		//cerr << curID << endl;
		SEEK = Targets.find(curID);
		if (SEEK == Targets.end()) {//no hit
			if (RemainderStrPos != -1) {
				MD->writeSelectiveStream(match[0], 0, RemainderStrPos);
				if (pairs > 1) {
					MD->writeSelectiveStream(match[1], 1, RemainderStrPos);
				}
			} else {// nothing written, but still need to delete
//				delete match; if (pairs > 1) {delete match2;}
			}
			continue;
		} 
		//serious work 
		//cerr << "H";
		idx = SEEK->second;
		if (b_doHD && newHD[idx] != "" && pairs>1) {
            match[0]->setHeader(newHD[idx] + "#1:0");
            match[1]->setHeader(newHD[idx] + "#2:0");
		}
		else if (b_doHD && newHD[idx] != "") {
			match[0]->setHeader(newHD[idx]); // single-end subsets were never renamed
		}
		MD->writeSelectiveStream(match[0], 0, outFilesIdx[idx]);
		if (pairs > 1) {
			MD->writeSelectiveStream(match[1], 1, outFilesIdx[idx]);
		}

		//finished, clean up to reduce search space (faster??)
		Targets.erase(SEEK);
		if (Targets.size() == 0 && RemainderStrPos==-1) {
			return;
		}
	}
	cerr << Targets.size() << " seqs remaining (not found in current file)\n";
}




//*******************************************
//*        DEREPLICATE OBJECT
//*******************************************

Dereplicate::ThreadLocalBatch& Dereplicate::get_or_create_thread_batch() {
 // The registry owns pending work even after a producer exits. The weak
 // cache detects destruction/address reuse without a reference-count RMW per read.
 thread_local std::weak_ptr<ThreadLocalBatch> cache;
 thread_local ThreadLocalBatch* cachedBatch = nullptr;
 if (cache.expired() || cachedBatch->owner != this) {
  auto batch = std::make_shared<ThreadLocalBatch>();
  batch->owner = this;
  batch->pending.reserve(kThreadBatchFlushThreshold);
  std::lock_guard<std::mutex> guard(batch_registry_mtx_);
  batch_registry_.push_back(batch);
  cache = batch;
  cachedBatch = batch.get();
 }
 return *cachedBatch;
}

size_t Dereplicate::tracker_size_total() const {
    size_t size = 0;
    for (const auto& shard : shards_) { size += shard->coarse ? (shard->coarse->outputsExact() ? shard->coarse->exactGroups().size() : shard->coarse->groups().size()) : shard->tracker.size(); }
    return size;
}

bool Dereplicate::apply_pending_entry_locked(const PendingDerepEntry& pending, HashDNA& tracker, bool deferring) {
    auto it = tracker.find(pending.srchSeq);
    if (it == tracker.end() && !pending.pass && !deferring) { return false; }
    bool fresh = it == tracker.end();
    if (fresh) { it = tracker.try_emplace(pending.srchSeq, std::make_unique<DNAuniqSet>()).first; }
    auto& group = it->second;
    if (coarse_ && deferring && !pending.pass &&
            (!coarse_->outputsExact() || coarse_->storesQualities() || coarse_->storesDiffs())) {
        // Reuse the ordinary deferred-prefix table, retaining the observations
        // needed for variant export as well as their counts. These reads never
        // supply the output representative; outputGroups attaches their evidence
        // only if a passing exact/prefix anchor exists.
        auto representative = group->find(pending.mergePos);
        if (representative == group->end()) {
            group->addNewDNAuniq(pending.dna, coarse_->storesQualities() ? pending.dna2 : nullptr,
                b_merge_pairs_derep_ ? pending.dna_merged : nullptr, pending.mergePos, -1, true, false);
            if (coarse_->storesQualities() || coarse_->storesDiffs()) {
                const auto& source = searchWithMerg && pending.dna_merged ? pending.dna_merged : pending.dna;
                auto reference = std::make_shared<const string>(coarse_->storesQualities() ? std::string_view(source->getSequence()) : source->getSeqPseudoView());
                auto mateReference = coarse_->storesQualities() && pending.dna2 ?
                    std::make_shared<const string>(pending.dna2->getSequence()) : nullptr;
                (*group)[pending.mergePos]->attachDiffDNA(std::make_shared<diffDNA>(reference, false,
                    coarse_->storesQualities(), mateReference));
            }
        }
        auto& dna = (*group)[pending.mergePos];
        if (coarse_->storesQualities() && pending.dna2) { dna->recordDiffPairQuality(*pending.dna, *pending.dna2, pending.sample_id); }
        else if (coarse_->storesQualities()) { dna->recordDiffQuality(*pending.dna, pending.sample_id); }
        else {
            const auto& source = searchWithMerg && pending.dna_merged ? pending.dna_merged : pending.dna;
            dna->recordDiffDNA(source->getSeqPseudoView(), pending.sample_id);
        }
        return true;
    }
    // Deferred failures are buffered separately. They become evidence only
    // if a passing anchor exists at output time; buffering is not acceptance.
    const bool evidenceOnly = deferring && !pending.pass;
    // Deferred failures can never supply an output representative. Keep only
    // R1 evidence and, when requested for output, merged quality evidence.
    // Their mates and unused merged products need not be retained in memory.
    group->addObservation(pending.dna, evidenceOnly ? nullptr : pending.dna2,
        evidenceOnly && !b_merge_pairs_derep_ ? nullptr : pending.dna_merged,
        pending.mergePos, pending.sample_id, b_derep_as_fasta_, b_merge_pairs_derep_, !evidenceOnly);
    return true;
}

void Dereplicate::apply_pending_entry(const PendingDerepEntry& pending, Shard& shard,
    const CoarseDereplication::Probe* probe, std::vector<CoarseDereplication::Update>* deferred) {
    // Called with a shard lock or exclusive work-pool ownership. Statistics are
    // finalized here, not while filtering still has this observation queued.
    auto& stats = shard.admission;
    stats.mergedKeys += searchWithMerg && bool(pending.dna_merged);
    // A non-passing read whose key already has a passing dereplicate is counted
    // right away. Otherwise a set with withheld mid-quality output is written
    // instead of counted (a pipeline back-maps it; never both), while a set
    // without such output is deferred: attached at output to the passing group
    // of its key (prefix mode: every non-passing read), or dropped, never
    // reaching .rest or the map on its own.
    if (!pending.pass && defersNonPassing()) {
        const bool known = exact_key_known(pending, shard);
        if (pending.withheld) {
            if (!known) { publish_uncounted(pending); ++stats.published; return; }
        } else if (prefixMode_) {
            if (deferCountsOnly_) { shard.rejectedTallies[pending.srchSeq].add(pending.sample_id); }
            else { apply_pending_entry_locked(pending, shard.rejected, true); }
            ++stats.deferred;
            return;
        } else if (!known) {
            bool kept;
            if (retainsVariants()) {
                // The reads themselves (already allocated), within a byte budget.
                auto bytes = [](const shared_ptr<DNA>& d) -> uint64_t {
                    return d ? sizeof(DNA) + 2 * uint64_t(d->mem_length()) + 2 * d->getId().size() : 0;
                };
                const uint64_t size = bytes(pending.dna) + bytes(pending.dna2) + bytes(pending.dna_merged);
                kept = deferredEvidenceBytes_.fetch_add(size, std::memory_order_relaxed) + size <= kDeferredEvidenceBytes;
                if (kept) {
                    shard.deferredEvidence[deferredKeyHash(pending.srchSeq)].push_back(
                        DeferredEvidence{pending.dna, pending.dna2, pending.dna_merged, pending.sample_id});
                }
                else { deferredEvidenceBytes_.fetch_sub(size, std::memory_order_relaxed); }
            }
            else {
                // Sample counts only, under a key hash, in a capped table.
                kept = shard.deferredCounts.add(deferredKeyHash(pending.srchSeq), pending.sample_id, deferredKeyCapPerShard_);
            }
            if (kept) { ++stats.deferred; }
            else { ++stats.unmatched; ++stats.deferredDropped; } // as 3.60: not counted
            return;
        }
    }
    bool admitted;
    if (shard.coarse && (!pending.srchSeq.empty() || !shard.coarse->storageOnly())) {
        const auto& source = searchWithMerg && pending.dna_merged ? pending.dna_merged : pending.dna;
        admitted = shard.coarse->add(pending.srchSeq, source->getSeqPseudoView(), pending.dna, pending.dna2,
            pending.dna_merged, pending.sample_id, pending.pass, searchLength,
            b_derep_as_fasta_, b_merge_pairs_derep_, pending.mergePos, probe, deferred) != CoarseDereplication::npos;
    } else {
        admitted = apply_pending_entry_locked(pending, shard.tracker, false);
    }
    if (!admitted) {
        if (pending.withheld) { publish_uncounted(pending); ++stats.published; }
        else { ++stats.unmatched; }
    }
    else if (pending.pass) { ++stats.high; }
    else if (pending.dna->isYellowQual()) { ++stats.mid; }
    else { ++stats.failed; }
}

void Dereplicate::publish_uncounted(const PendingDerepEntry& pending) {
    if (!publisher_) { throw std::logic_error("Dereplicate: a withheld record set has no publisher"); }
    publisher_(pending.dna, pending.dna2, pending.route);
}

void Dereplicate::flushPending() {
    LifecyclePause pause(*this);
    flush_all_batches();
}

bool Dereplicate::DeferredCounts::add(uint64_t key, int sample, size_t cap) {
    if (sample < 0) { return true; } // no sample: nothing to count (diffDNA ignores it as well)
    auto bump = [sample](std::vector<std::pair<int, uint64_t>>& counts) {
        for (auto& entry : counts) { if (entry.first == sample) { ++entry.second; return; } }
        counts.emplace_back(sample, 1);
    };
    if (auto several = multi.find(key); several != multi.end()) { bump(several->second); return true; }
    auto one = single.find(key);
    if (one == single.end()) {
        if (size() >= cap) { return false; }
        single.emplace(key, Single{sample, 1});
        return true;
    }
    if (one->second.sample == sample && one->second.copies < std::numeric_limits<uint32_t>::max()) {
        ++one->second.copies;
        return true;
    }
    std::vector<std::pair<int, uint64_t>> counts{{one->second.sample, one->second.copies}};
    single.erase(one);
    bump(counts);
    multi.emplace(key, std::move(counts));
    return true;
}

bool Dereplicate::exact_key_known(const PendingDerepEntry& pending, const Shard& shard) const {
    if (shard.coarse && (!pending.srchSeq.empty() || !shard.coarse->storageOnly())) {
        return shard.coarse->hasExactKey(pending.srchSeq);
    }
    return shard.tracker.find(pending.srchSeq) != shard.tracker.end();
}

void Dereplicate::flush_batch(ThreadLocalBatch& batch) {
	if (batch.pending.empty()) { return; }
	uint64_t searchBases = 0;
	for (const auto& entry : batch.pending) { searchBases += entry.srchSeq.size(); }
	srch_seq_total_len_.fetch_add(searchBases, std::memory_order_relaxed);
	srch_seq_count_.fetch_add(batch.pending.size(), std::memory_order_relaxed);
	struct ClearBatch { ThreadLocalBatch& batch; ~ClearBatch() { batch.pending.clear(); } } clear{batch};
    if (parallelGlobalCoarse_) {
        flush_global_coarse(batch);
        return;
    }
	thread_local std::vector<std::vector<size_t>> grouped;
	grouped.resize(shards_.size());
	for (auto& indices : grouped) { indices.clear(); }
	for (size_t i = 0; i < batch.pending.size(); ++i) { grouped[batch.pending[i].shard].push_back(i); }
    if (coarseWorkPool_) {
        std::vector<PartitionWorkPool::Work> work;
        work.reserve(shards_.size());
        for (size_t i = 0; i < grouped.size(); ++i) {
            if (!grouped[i].empty()) { work.push_back({i, grouped[i].size()}); }
        }
        // Capture the submitting thread's grouping explicitly. A helper must
        // not resolve the thread_local variable to its own current batch.
        auto* indices = &grouped;
        coarseWorkPool_->run(work, [this, &batch, &work, indices](size_t task) {
            auto& shard = *shards_[work[task].partition];
            // Pool ownership serializes this index, without a shard mutex wait.
            for (size_t index : (*indices)[work[task].partition]) {
                apply_pending_entry(batch.pending[index], shard);
            }
        });
        return;
    }
	// Shards are independent: apply to every shard that is free right away and
	// only then wait for the busy ones, instead of queueing behind other
	// workers in a fixed shard order.
	thread_local std::vector<size_t> busy;
	busy.clear();
	for (size_t i = 0; i < shards_.size(); ++i) {
		if (grouped[i].empty()) { continue; }
		auto& shard = *shards_[i];
		std::unique_lock<std::mutex> guard(shard.mutex, std::try_to_lock);
		if (!guard.owns_lock()) { busy.push_back(i); continue; }
        for (size_t index : grouped[i]) {
            apply_pending_entry(batch.pending[index], shard);
        }
	}
	for (size_t i : busy) {
		auto& shard = *shards_[i];
		std::lock_guard<std::mutex> guard(shard.mutex);
        for (size_t index : grouped[i]) {
            apply_pending_entry(batch.pending[index], shard);
        }
	}
}

void Dereplicate::flush_global_coarse(ThreadLocalBatch& batch) {
    auto& shard = *shards_.front();
    auto& coarse = *shard.coarse;
    thread_local CoarseDereplication::SearchScratch scratch;
    thread_local std::vector<CoarseDereplication::Probe> probes;
    thread_local std::vector<CoarseDereplication::Update> updates;
    probes.assign(batch.pending.size(), CoarseDereplication::Probe{});
    updates.clear();
    struct ClearUpdates { ~ClearUpdates() { updates.clear(); } } clearUpdates;
    {
        // Phase 1: candidate searches (most of the work) against a snapshot.
        std::shared_lock<std::shared_mutex> snapshot(shard.snapshot);
        for (size_t i = 0; i < batch.pending.size(); ++i) {
            const auto& pending = batch.pending[i];
            if (prefixMode_ && !pending.pass) { continue; } // deferred prefix evidence, no search
            probes[i] = coarse.probe(pending.srchSeq, pending.dna, pending.dna2, pending.pass, searchLength, scratch);
        }
    }
    {
        // Phase 2: commit in batch order. Exact keys are re-checked and only
        // anchors created after each snapshot are searched again, so the
        // result equals committing this batch serially at this point.
        std::lock_guard<std::mutex> guard(shard.mutex);
        std::unique_lock<std::shared_mutex> exclusive(shard.snapshot);
        for (size_t i = 0; i < batch.pending.size(); ++i) {
            apply_pending_entry(batch.pending[i], shard, &probes[i], &updates);
        }
    }
    // Phase 3: representative, count and variant updates touch only their own
    // target; updates of one target stay in batch order under its stripe.
    for (const auto& update : updates) {
        const size_t stripe = (reinterpret_cast<uintptr_t>(update.target) >> 6) % coarseUpdateStripes_.size();
        std::lock_guard<std::mutex> guard(coarseUpdateStripes_[stripe]);
        coarse.apply(update);
    }
}

void Dereplicate::flush_all_batches() {
	drain_sealed_blocks();
	std::vector<std::shared_ptr<ThreadLocalBatch>> batches;
	{ std::lock_guard<std::mutex> guard(batch_registry_mtx_); batches = batch_registry_; }
	for (auto& batch : batches) {
		flush_batch(*batch);
		free_returned_reads(*batch);
	}
}

Dereplicate::Dereplicate(OptContainer* cmdArgs, Filters* mf):
        barcode_number_to_sample_id_(0), b_usearch_fmt(true), b_singleLine(true), b_pairedInput(false),
        minCopyRules(1, DerepCopyRule{0, 0}), minCopiesStr("0"), //default accepts every derep
		minCopiesFastFailThreshold(-1),
		tmpCnt(0), curBCoffset(0), b_derep_as_fasta_(true), b_derepPerSR(false),
		b_wroteMapHD(false), b_merge_pairs_derep_(false),merger(nullptr),
		searchLength(150),
		mapF(""), outHQf(""), outHQf_p2(""), outRest(""),
		mainFilter(mf), searchWithMerg(true)
{
    const uint64_t identity = CoarseDereplication::parseIdentity(cmdArgs->get("-derepIdentity", "100"));
    const unsigned homopolymer = CoarseDereplication::parseHomopolymer(cmdArgs->get("-derepSeedHomopolymer", "4"));
    const string diffOption = cmdArgs->get("-derepStoreDiffs", "0");
    if (diffOption != "0" && diffOption != "1") { throw std::invalid_argument("-derepStoreDiffs must be 0 or 1"); }
    const string subclusterOption = cmdArgs->get("-derepSubclusterFasta", "0");
    if (subclusterOption != "0" && subclusterOption != "1") {
        throw std::invalid_argument("-derepSubclusterFasta must be 0 or 1");
    }
    writeSubclusters_ = subclusterOption == "1";
    if (writeSubclusters_ && diffOption != "1") {
        throw std::invalid_argument("-derepSubclusterFasta 1 requires -derepStoreDiffs 1");
    }
    const string reassignOption = cmdArgs->get("-derepReassign", "0");
    if (reassignOption != "0" && reassignOption != "1") { throw std::invalid_argument("-derepReassign must be 0 or 1"); }
    const string clusterOption = cmdArgs->get("-derepCoarseClusters", "0");
    if (clusterOption != "0" && clusterOption != "1") { throw std::invalid_argument("-derepCoarseClusters must be 0 or 1"); }
    // Explicit final reassignment also requests coarse cluster output.
    const bool exactOutput = clusterOption != "1" && reassignOption != "1";
    const string qualityOption = cmdArgs->get("-derepStoreQuals", "0");
    if (qualityOption != "0" && qualityOption != "1") { throw std::invalid_argument("-derepStoreQuals must be 0 or 1"); }
    if (qualityOption == "1") {
        if (!cmdArgs->get("-i_fna").empty() && cmdArgs->get("-i_qual").empty()) {
            throw std::invalid_argument("-derepStoreQuals requires complete input qualities; FASTA alone is insufficient");
        }
        const auto fastaFiles = mf->getFastaFiles(), qualityFiles = mf->getQualFiles();
        for (size_t i = 0; i < fastaFiles.size(); ++i) {
            if (!fastaFiles[i].empty() && (i >= qualityFiles.size() || qualityFiles[i].empty())) {
                throw std::invalid_argument("-derepStoreQuals requires a quality file for every FASTA input");
            }
        }
    }
    if (qualityOption == "1") {
        b_pairedInput = mf->isPaired() > 1;
        // Exact pair variants keep the two processed mates and their own scores.
        searchWithMerg = false;
        for (const auto* option : {"-merge_pairs_derep", "-merge_pairs_filter", "-merge_pairs_demulti"}) {
            if (cmdArgs->get(option) == "1") { throw std::invalid_argument(std::string("-derepStoreQuals requires unmerged input; disable ") + option); }
        }
    }
    const string globalOption = cmdArgs->get("-derepGlobal", "0");
    if (globalOption != "0" && globalOption != "1") { throw std::invalid_argument("-derepGlobal must be 0 or 1"); }
    const bool coarseMode = identity < CoarseDereplication::identityScale || diffOption == "1" || reassignOption == "1" || qualityOption == "1";
    const string prefixOption = cmdArgs->contains("-derepPrefix") ? (*cmdArgs)["-derepPrefix"] : "auto";
    if (prefixOption != "auto" && prefixOption != "0" && prefixOption != "1") {
        throw std::invalid_argument("-derepPrefix must be auto, 0, or 1");
    }
	const OptContainer::RuntimeOptions& runtimeOpts = cmdArgs->runtime();
	outfile = (*cmdArgs)["-o_dereplicate"];
	if (runtimeOpts.noSearchWithMerge) {
		searchWithMerg = false;
	}

	string baseOF = std::filesystem::path(outfile).replace_extension().string();
    if (coarseMode && diffOption == "1" && outfile == baseOF + ".diff") {
        throw std::invalid_argument("The dereplicate output cannot be its own .diff sidecar; choose a FASTA/FASTQ output filename");
    }
	mapF = baseOF + ".map";
	outHQf = baseOF + ".1.hq.fq";
	outHQf_p2 = baseOF + ".2.hq.fq";
	outRest = outfile + ".rest";
    //clean up eventually existing files
	remove(mapF.c_str());	remove(outHQf.c_str());
	remove(outHQf_p2.c_str());	remove(outRest.c_str());

	// Reserve hash map capacity to avoid costly rehashing during heavy dereplication.
	// Allow user to override via -derep_reserve; default to 500k slots.
	size_t reserveSize = coarseMode ? 0 : runtimeOpts.derepReserve;
// Independent sequence keys can update concurrently. Keep each key in one
	// shard so occurrence counts, quality evidence, and selection stay atomic.
    // Approximate groups may stay independent. Exact-prefix and explicitly global
    // clustering/reassignment retain one index so cross-key matching stays global.
    // Explicit coarse-cluster output also needs the global index: approximate
    // matches in different partitions would never be grouped, and the
    // partitioned mode exists to save memory for exact-parent output.
    const bool globalCoarse = runtimeOpts.threads <= 1 || globalOption == "1" ||
        reassignOption == "1" || clusterOption == "1" || identity == CoarseDereplication::identityScale;
    const size_t shardCount = coarseMode && globalCoarse ? 1 : 64;
	shards_.reserve(shardCount);
	for (size_t i = 0; i < shardCount; ++i) {
		auto shard = std::make_unique<Shard>();
		shard->tracker.reserve(reserveSize / shardCount + (reserveSize % shardCount != 0));
        if (coarseMode) {
            shard->coarse = std::make_unique<CoarseDereplication>(identity, homopolymer, diffOption == "1", reassignOption == "1", qualityOption == "1", exactOutput);
        }
		shards_.push_back(std::move(shard));
	}
    if (coarseMode) { coarse_ = shards_.front()->coarse.get(); }
    if (coarseMode && shardCount > 1) { coarseWorkPool_ = std::make_unique<PartitionWorkPool>(shardCount); }
    parallelGlobalCoarse_ = coarseMode && shardCount == 1 && runtimeOpts.threads > 1;
    orderedAdmission_ = runtimeOpts.threads > 1;
    if (orderedAdmission_ && !coarseWorkPool_ && !parallelGlobalCoarse_) { orderedPool_ = std::make_unique<PartitionWorkPool>(shardCount); }
    // Enough queued blocks to keep the committer busy while every other worker
    // processes its next block; queued blocks keep their reads in memory.
    commitBacklog_ = std::max<size_t>(4, 2 * size_t(runtimeOpts.threads));
    // A few blocks' worth: enough for every lane to have work, small enough
    // that the committer's final drain stays short.
    updateBacklog_ = 256 * size_t(runtimeOpts.threads);


	if (runtimeOpts.hasDereSizeFmt && runtimeOpts.disableUsearchSizeFmt) {
		b_usearch_fmt = false;
	}
	if (runtimeOpts.derepPerSR) {
		b_derepPerSR = true;
	}

	if (runtimeOpts.hasMinDerepCopies) {
		minCopiesStr = runtimeOpts.minDerepCopies;
		string x = runtimeOpts.minDerepCopies;
		std::replace(x.begin(), x.end(), ';', ','); // "8:1;3:4" is accepted like "8:1,3:4"
		minCopyRules.clear();

		vector<string> xs = splitByCommas(x, ',');
		for (size_t i = 0; i < xs.size(); i++) {
			if (xs[i].empty()) { continue; }
			vector<string> tmp = splitByComma(xs[i], false, ':');
			if (tmp.size() > 2) {
				cerr << "Derep string " << xs[i] << " has to be two integers seqparated by \":\"\n"; exit(623);
			}
			// A bare number counts copies across all samples; "X:Y" needs X copies in each of Y samples.
			DerepCopyRule rule{ atoll(tmp[0].c_str()), 0 };
			if (tmp.size() == 2) {
				rule.samples = atoi(tmp[1].c_str());//number of samples
				if (rule.samples <= 0) { cerr << "wrong derplicate position (<=0):" << xs[i] << endl; exit(312); }
				if (rule.samples > 3000) { cerr << "too large derplicate position (>3000):" << xs[i] << endl; exit(313); }
			}
			minCopyRules.push_back(rule);
		}
		if (minCopyRules.empty()) { minCopyRules.push_back(DerepCopyRule{ 0, 0 }); }
	}

	// Smallest total count any rule can pass with (-1: some rule accepts everything).
	minCopiesFastFailThreshold = -1;
	for (const DerepCopyRule& rule : minCopyRules) {
		const int64_t needed = rule.copies * std::max(1, rule.samples);
		if (needed <= 0) {
			minCopiesFastFailThreshold = -1;
			break;
		}
		if (minCopiesFastFailThreshold == -1 || needed < minCopiesFastFailThreshold) {
			minCopiesFastFailThreshold = needed;
		}
	}
	// set up output format of dereplication:
	if (runtimeOpts.hasDerepFormat) {
		auto res = runtimeOpts.derepFormat;
		if (strcmp(res.c_str(), "fa") == 0) {
			b_derep_as_fasta_ = true;
		}
		else if (strcmp(res.c_str(), "fq") == 0) {
			b_derep_as_fasta_ = false;
		}
		else {
			std::cerr << "Invalid input for option -derep_format. If specified, -derep_format can take either \"fa\" for derep_as_fasta_ and \"fq\" for fastq. ";
			exit(77);
		}
	}

	// Flag for merging paired end fsatq reads
	if (runtimeOpts.hasMergePairsDerep && runtimeOpts.mergePairsDerep) {
		b_merge_pairs_derep_ = true;
	}

	if (runtimeOpts.hasDerepSrchLen) {
		const string& searchLenArg = (*cmdArgs)["-derepSrchLen"];
		size_t parsedChars = 0;
		long long parsedLen = 0;
		try {
			parsedLen = std::stoll(searchLenArg, &parsedChars);
		}
		catch (...) {
			cerr << "Invalid search length specified for -derepSrchLen. Search length must be -1 or a positive integer. Exiting.\n";
			exit(78);
		}
		if (parsedChars != searchLenArg.size() || parsedLen == 0 || parsedLen < -1 || parsedLen > std::numeric_limits<int>::max()) {
			cerr << "Invalid search length specified for -derepSrchLen. Search length must be -1 or a positive integer. Exiting.\n";
			exit(78);
		}
		searchLength = static_cast<int>(parsedLen);
	}
    const int fixedLength = mf->fixedDereplicationLength(searchWithMerg);
    const int clusteringCap = mf->dereplicationLengthLimit();
    const bool explicitCap = mf->hasExplicitDerepSearchLength() && searchLength > 0;
    prefixMode_ = (!coarse_ || coarse_->outputsExact()) && (prefixOption == "1" || (prefixOption == "auto" && !explicitCap && fixedLength == 0));
    if (prefixMode_) {
        searchLength = clusteringCap;
    } else if (coarse_) {
        if (!explicitCap) { searchLength = -1; }
        if (clusteringCap > 0 && (searchLength < 0 || searchLength > clusteringCap)) { searchLength = clusteringCap; }
    } else if (prefixOption == "auto") {
        if (!explicitCap && fixedLength > 0) { searchLength = fixedLength; }
        if (clusteringCap > 0 && (searchLength < 0 || searchLength > clusteringCap)) {
            searchLength = clusteringCap;
        }
    }
    cerr << "Dereplication matching: " << (coarse_ ? cmdArgs->get("-derepIdentity", "100") + "% streaming substitution clusters" : (prefixMode_ ? "100% prefix" : "exact length"))
         << (searchLength > 0 ? " (search cap " + to_string(searchLength) + ")" : " (full logical sequence)") << "\n";
    if (coarse_) {
        cerr << "Coarse dereplication: " << shards_.size() << (shards_.size() > 1 ?
            " independent partitions; no final consolidation\n" : " global index\n");
        if (coarseWorkPool_) { cerr << "Coarse scheduling: shared partition work pool (existing workers)\n"; }
        if (parallelGlobalCoarse_) { cerr << "Coarse scheduling: concurrent snapshot searches, ordered commits\n"; }
    }
    // Prefix mode: FASTA output collects no quality evidence from deferred
    // reads, so only their sample counts need to be kept; retained variants and
    // FASTQ output keep the read (as before). Exact-key mode keeps counts only.
    deferCountsOnly_ = b_derep_as_fasta_ && !(coarse_ && (coarse_->storesQualities() || coarse_->storesDiffs()));
    deferredKeyCapPerShard_ = std::max<size_t>(1, kDeferredKeyCap / std::max<size_t>(1, shards_.size()));
}

Dereplicate::~Dereplicate() {
	if (merger != nullptr) { delete merger; merger = nullptr; }
}


Dereplicate::ActiveAddGuard::ActiveAddGuard(Dereplicate& owner, ThreadLocalBatch& batch)
	: owner(owner), batch(batch) {
	for (;;) {
		if (owner.lifecycle_pause_.load()) {
			std::unique_lock<std::mutex> lock(owner.lifecycle_wait_mtx_);
			owner.lifecycle_wait_cv_.wait(lock, [&owner] { return !owner.lifecycle_pause_.load(); });
		}
		// Sequential consistency pairs this announcement/recheck with the
		// snapshot's pause/store and activity scan. Either it sees this producer,
		// or this producer sees the pause before touching any pending entries.
		batch.active.store(true);
		if (!owner.lifecycle_pause_.load()) { return; }
		owner.release_batch_activity(batch);
	}
}

Dereplicate::ActiveAddGuard::~ActiveAddGuard() { owner.release_batch_activity(batch); }

void Dereplicate::release_batch_activity(ThreadLocalBatch& batch) {
	batch.active.store(false);
	if (lifecycle_pause_.load()) {
		// Synchronize notification with the waiter's predicate check so the
		// final producer cannot notify just before the waiter goes to sleep.
		std::lock_guard<std::mutex> lock(lifecycle_wait_mtx_);
		lifecycle_wait_cv_.notify_all();
	}
}

Dereplicate::LifecyclePause::LifecyclePause(Dereplicate& owner)
	: owner(owner), operationLock(owner.lifecycle_operation_mtx_) {
	std::unique_lock<std::mutex> lock(owner.lifecycle_wait_mtx_);
	std::vector<std::shared_ptr<ThreadLocalBatch>> batches;
	{
		std::lock_guard<std::mutex> registryLock(owner.batch_registry_mtx_);
		// Copy before publishing the pause: allocation failure leaves producers running.
		batches = owner.batch_registry_;
		owner.lifecycle_pause_.store(true);
	}
	// Newly registered producers see the pause before touching their batches.
	// Release the registry lock before waiting/reacquiring the condition mutex.
	owner.lifecycle_wait_cv_.wait(lock, [&batches] {
		return std::all_of(batches.begin(), batches.end(), [](const auto& batch) {
			return !batch->active.load();
		});
	});
}

Dereplicate::LifecyclePause::~LifecyclePause() {
	std::lock_guard<std::mutex> lock(owner.lifecycle_wait_mtx_);
	owner.lifecycle_pause_.store(false);
	owner.lifecycle_wait_cv_.notify_all();
}

bool Dereplicate::addDNA(shared_ptr<DNA> dna, shared_ptr<DNA> dna2, Filters* filter,
    const DerepTierRoute* route, bool* queued) {
    if (queued) { *queued = false; }
    if (coarse_ && coarse_->storesQualities() && bool(dna2) != b_pairedInput) {
        throw std::invalid_argument("-derepStoreQuals requires complete pairs for paired input, or consistently single-end input");
    }
	ThreadLocalBatch& threadBatch = get_or_create_thread_batch();
	ActiveAddGuard activeGuard(*this, threadBatch);
    if (coarse_ && coarse_->storesQualities() && dna2 && dna2->mem_length() == 0) {
        skippedEmptyMates_.fetch_add(1, std::memory_order_relaxed);
        return false;
    }
	//1st build hash of DNA
	if (!dna->getBarcodeDetected()) {
		return false;
	}
	// Get copy of sequence (might have already been modified)
	//string seq = dna->getSeqPseudo();

	int sample_id = dna->getBCnumber();
	bool pass = dna->isGreenQual();
	//deactivate this for now..

	int MrgPos1 = dna->merge_seed_pos_;
	//int MrgPos1 = -1;
	shared_ptr<DNA> dna_merged = nullptr;
	const string* srchSeqPtr = nullptr;
	thread_local string srchSeqBuffer;

	/*
	* findSeedForMerge requires own thread..
	if (searchWithMerg && merger != nullptr && dna->merge_seed_pos_ == -1) {
		MrgPos1 = merger->findSeedForMerge(dna, dna2);
	}*/

    // Only attempt to merge if a merger has been activated and pair is present
	if ((searchWithMerg || b_merge_pairs_derep_) && dna2 != nullptr && dna->merge_seed_pos_ != -1 && merger != nullptr) {
		//MrgPos1 = dna2->merge_seed_pos_;
		dna_merged = merger->merge(dna, dna2);
        Filters* mergedFilter = filter ? filter : mainFilter;
        if (dna_merged && !mergedFilter->checkMergedDNA(dna_merged)) {
            dna_merged.reset();
        }
        // Only an accepted consensus marks its mates as merged (length and
        // overlap errors feed representative selection and the report).
        if (dna_merged) { dna->adoptMergeStats(*dna_merged); dna2->adoptMergeStats(*dna_merged); }
	}

    if (dna_merged && searchWithMerg){
		const std::string_view mergedSeq = dna_merged->getSeqPseudoView();
		const size_t searchLen = (searchLength == -1)
			? mergedSeq.size()
			: std::min(static_cast<size_t>(searchLength), mergedSeq.size());
		srchSeqBuffer.assign(mergedSeq.data(), searchLen);
		srchSeqPtr = &srchSeqBuffer;
	} else {
		const std::string_view seq = dna->getSeqPseudoView();
		const size_t searchLen = (searchLength == -1)
			? seq.size()
			: std::min(static_cast<size_t>(searchLength), seq.size());
		srchSeqBuffer.assign(seq.data(), searchLen);
		srchSeqPtr = &srchSeqBuffer;
	}
	const string& srchSeq = *srchSeqPtr;
	if ((prefixMode_ || (coarse_ && !coarse_->storageOnly())) && srchSeq.empty()) {
        skippedEmptyKeys_.fetch_add(1, std::memory_order_relaxed);
        return false;
    }

    // Exact dereplication is keyed by R1/search sequence, including in coarse
    // mode. Different R2s must reach the same owner to compete as best seeds.
    const size_t keyHash = coarse_ && coarse_->storageOnly() && shards_.size() > 1 ?
        coarse_->storagePartitionHash(srchSeq) : robin_hood::hash<string>{}(srchSeq);
	threadBatch.pending.push_back(PendingDerepEntry{
		srchSeq,
		dna,
		dna2,
		dna_merged,
		MrgPos1,
		sample_id,
		pass,
		keyHash % shards_.size()
		});
	if (route) { threadBatch.pending.back().withheld = true; threadBatch.pending.back().route = *route; }
	if (queued) { *queued = true; }
	// Inside an ordered block, entries are applied when the block is sealed.
	bool flushNow = threadBatch.pending.size() >= kThreadBatchFlushThreshold;

	if (threadBatch.block == kNoBlock && (flushNow || lifecycle_pause_.load(std::memory_order_acquire))) {
		flush_batch(threadBatch);
	}

	return pass;
}

uint64_t Dereplicate::reserveBlock() {
	return orderedAdmission_ ? nextBlock_.fetch_add(1, std::memory_order_relaxed) : kNoBlock;
}

void Dereplicate::beginBlock(uint64_t block) {
	if (block == kNoBlock) { return; }
	ThreadLocalBatch& batch = get_or_create_thread_batch();
	if (!batch.pending.empty()) {
		// Entries added outside any block keep their unordered path.
		ActiveAddGuard activeGuard(*this, batch);
		flush_batch(batch);
	}
	batch.block = block;
}

void Dereplicate::account_search_keys(const std::vector<PendingDerepEntry>& entries) {
	uint64_t searchBases = 0;
	for (const auto& entry : entries) { searchBases += entry.srchSeq.size(); }
	srch_seq_total_len_.fetch_add(searchBases, std::memory_order_relaxed);
	srch_seq_count_.fetch_add(entries.size(), std::memory_order_relaxed);
}

void Dereplicate::endBlock() {
	ThreadLocalBatch& batch = get_or_create_thread_batch();
	if (batch.block == kNoBlock) { return; }
	const uint64_t block = batch.block;
	batch.block = kNoBlock;
	free_returned_reads(batch);
	bool queued = false;
	{
		ActiveAddGuard activeGuard(*this, batch);
		std::vector<PendingDerepEntry> entries;
		entries.swap(batch.pending);
		batch.pending.reserve(kThreadBatchFlushThreshold);
		account_search_keys(entries);
		if (parallelGlobalCoarse_ && !entries.empty()) {
			// The candidate search (most of the work) runs here, concurrently with
			// other workers, against a snapshot; the ordered commit reuses it.
			auto& shard = *shards_.front();
			thread_local CoarseDereplication::SearchScratch scratch;
			for (size_t begin = 0; begin < entries.size(); begin += kProbeChunk) {
				// The shared lock prefers readers: step aside for a waiting commit.
				while (commitWaiting_.load(std::memory_order_acquire)) { commitWaiting_.wait(true, std::memory_order_acquire); }
				std::shared_lock<std::shared_mutex> snapshot(shard.snapshot);
				const size_t end = std::min(entries.size(), begin + kProbeChunk);
				for (size_t i = begin; i < end; ++i) {
					auto& pending = entries[i];
					if (prefixMode_ && !pending.pass) { continue; } // deferred prefix evidence, no search
					pending.probe = shard.coarse->probe(pending.srchSeq, pending.dna, pending.dna2, pending.pass, searchLength, scratch);
				}
			}
		}
		std::vector<SealedBlock> ready;
		std::unique_lock<std::mutex> submitLock(submitMutex_, std::defer_lock);
		bool commit = false;
		{
			std::lock_guard<std::mutex> lock(sequenceMutex_);
			sealed_.emplace(block, SealedBlock{std::move(entries), &batch});
			while (!sealed_.empty() && sealed_.begin()->first == nextSubmit_) {
				ready.push_back(std::move(sealed_.begin()->second));
				sealed_.erase(sealed_.begin());
				++nextSubmit_;
			}
			// An earlier block is still being processed: its worker (or whoever
			// completes the gap) applies this block after it.
			if (ready.empty()) { return; }
			// Queued, or the submit lock taken, before releasing sequenceMutex_,
			// so whoever takes the next blocks can only enqueue them after these.
			if (parallelGlobalCoarse_) { commit = queue_global_commits(ready); }
			else { submitLock.lock(); }
		}
		if (!parallelGlobalCoarse_) { apply_blocks_ordered(ready, submitLock, &batch); }
		else if (commit) { commit_global_blocks(); }
		else { queued = true; }
	}
	// Outside the activity guard: a lifecycle pause need not wait for this.
	if (queued) { wait_for_commit_backlog(); }
}

bool Dereplicate::queue_global_commits(std::vector<SealedBlock>& blocks) {
	std::lock_guard<std::mutex> lock(commitMutex_);
	if (commitError_) { std::rethrow_exception(commitError_); }
	for (auto& block : blocks) { commitQueue_.push_back(std::move(block)); }
	if (committing_) { return false; }
	committing_ = true;
	return true;
}

void Dereplicate::commit_global_blocks() {
	// Only one committer at a time (committing_). It drains the queue, including
	// blocks queued meanwhile, and every update lane before it returns, so
	// queued work always has an active committer (which a lifecycle pause
	// waits for) and helpers only ever run while it is active.
	std::unique_lock<std::mutex> lock(commitMutex_);
	while (!commitError_) {
		if (!commitQueue_.empty()) {
			SealedBlock block = std::move(commitQueue_.front());
			commitQueue_.pop_front();
			commitSpace_.notify_all();
			lock.unlock();
			try {
				commit_global_block(block);
				return_reads(block.owner, std::move(block.entries));
			} catch (...) {
				lock.lock();
				fail_commits_locked(std::current_exception());
				continue;
			}
			lock.lock();
		} else if (!apply_ready_lane(lock)) {
			if (!pendingUpdates_) { break; }
			commitSpace_.wait(lock); // the remaining lanes are being applied by helpers
		}
	}
	// After a failure, lanes still being applied must finish before returning.
	commitSpace_.wait(lock, [this] { return !pendingUpdates_; });
	committing_ = false;
	commitSpace_.notify_all();
	if (commitError_) { std::rethrow_exception(commitError_); }
}

void Dereplicate::commit_global_block(SealedBlock& block) {
	auto& shard = *shards_.front();
	thread_local std::vector<CoarseDereplication::Update> updates;
	updates.clear();
	{
		// Commit in block order. Exact keys are re-checked and only anchors
		// created after each probe's snapshot are searched again, so the
		// result equals applying these entries serially at this point.
		std::lock_guard<std::mutex> guard(shard.mutex);
		commitWaiting_.store(true, std::memory_order_release);
		std::unique_lock<std::shared_mutex> exclusive(shard.snapshot);
		commitWaiting_.store(false, std::memory_order_release);
		commitWaiting_.notify_all();
		for (auto& pending : block.entries) { apply_pending_entry(pending, shard, &pending.probe, &updates); }
	}
	post_updates(updates, block.owner);
}

void Dereplicate::post_updates(std::vector<CoarseDereplication::Update>& updates, ThreadLocalBatch* owner) {
	if (updates.empty()) { return; }
	// Updates only touch their own target, and a target always maps to the
	// same lane: grouping keeps each target's updates in commit order.
	thread_local std::array<std::vector<CoarseDereplication::Update>, kUpdateLanes> byLane;
	for (auto& update : updates) {
		byLane[(reinterpret_cast<uintptr_t>(update.target) >> 6) % kUpdateLanes].push_back(std::move(update));
	}
	updates.clear();
	std::lock_guard<std::mutex> lock(commitMutex_);
	for (size_t id = 0; id < kUpdateLanes; ++id) {
		if (byLane[id].empty()) { continue; }
		auto& lane = updateLanes_[id];
		pendingUpdates_ += byLane[id].size();
		if (lane.queue.empty() && !lane.running) { readyLanes_.push_back(id); }
		lane.queue.push_back(UpdateGroup{std::move(byLane[id]), owner});
		byLane[id].clear();
	}
	commitSpace_.notify_all();
}

bool Dereplicate::apply_ready_lane(std::unique_lock<std::mutex>& lock) {
	if (readyLanes_.empty()) { return false; }
	const size_t id = readyLanes_.front();
	readyLanes_.pop_front();
	auto& lane = updateLanes_[id];
	lane.running = true;
	std::deque<UpdateGroup> work;
	work.swap(lane.queue);
	lock.unlock();
	size_t count = 0;
	std::exception_ptr error;
	try {
		auto& coarse = *shards_.front()->coarse;
		for (const auto& group : work) {
			for (const auto& update : group.updates) {
				// Stripes serialize with unordered flushes (flush_global_coarse).
				const size_t stripe = (reinterpret_cast<uintptr_t>(update.target) >> 6) % coarseUpdateStripes_.size();
				std::lock_guard<std::mutex> guard(coarseUpdateStripes_[stripe]);
				coarse.apply(update);
			}
		}
	} catch (...) { error = std::current_exception(); }
	for (auto& group : work) {
		count += group.updates.size();
		return_reads(group.owner, std::move(group.updates));
	}
	work.clear();
	lock.lock();
	lane.running = false;
	pendingUpdates_ -= count;
	if (error) { fail_commits_locked(error); }
	else if (!lane.queue.empty()) { readyLanes_.push_back(id); }
	commitSpace_.notify_all();
	return true;
}

void Dereplicate::fail_commits_locked(std::exception_ptr error) {
	if (!commitError_) { commitError_ = error; }
	commitQueue_.clear();
	for (auto& lane : updateLanes_) {
		for (const auto& group : lane.queue) { pendingUpdates_ -= group.updates.size(); }
		lane.queue.clear();
	}
	readyLanes_.clear();
	commitSpace_.notify_all();
}

void Dereplicate::wait_for_commit_backlog() {
	// Queued blocks and updates hold their reads: while either backlog is
	// exceeded, apply updates rather than start a new block.
	std::unique_lock<std::mutex> lock(commitMutex_);
	while (committing_ && !commitError_ &&
		(commitQueue_.size() > commitBacklog_ || pendingUpdates_ > updateBacklog_)) {
		if (!apply_ready_lane(lock)) { commitSpace_.wait(lock); }
	}
	if (commitError_) { std::rethrow_exception(commitError_); }
}

// Freeing a read in another thread's malloc arena contends with that thread's
// own allocations, so committed entries and applied updates go back to the
// producer that allocated them; it frees them before its next block.
void Dereplicate::return_reads(ThreadLocalBatch* owner, std::vector<PendingDerepEntry>&& entries) {
	if (!owner || entries.empty()) { return; }
	std::lock_guard<std::mutex> lock(owner->returnedMutex);
	owner->returnedEntries.push_back(std::move(entries));
}

void Dereplicate::return_reads(ThreadLocalBatch* owner, std::vector<CoarseDereplication::Update>&& updates) {
	if (!owner || updates.empty()) { return; }
	std::lock_guard<std::mutex> lock(owner->returnedMutex);
	owner->returnedUpdates.push_back(std::move(updates));
}

void Dereplicate::free_returned_reads(ThreadLocalBatch& batch) {
	std::vector<std::vector<PendingDerepEntry>> entries;
	std::vector<std::vector<CoarseDereplication::Update>> updates;
	{
		std::lock_guard<std::mutex> lock(batch.returnedMutex);
		entries.swap(batch.returnedEntries);
		updates.swap(batch.returnedUpdates);
	}
}

void Dereplicate::apply_blocks_ordered(std::vector<SealedBlock>& blocks,
	std::unique_lock<std::mutex>& submitLock, ThreadLocalBatch* self) {
	struct Task { size_t block, shard; std::vector<size_t> indices; };
	std::vector<Task> tasks;
	std::vector<PartitionWorkPool::Work> work;
	std::vector<std::vector<size_t>> grouped(shards_.size());
	for (size_t b = 0; b < blocks.size(); ++b) {
		for (size_t i = 0; i < blocks[b].entries.size(); ++i) { grouped[blocks[b].entries[i].shard].push_back(i); }
		for (size_t s = 0; s < grouped.size(); ++s) {
			if (grouped[s].empty()) { continue; }
			work.push_back({s, grouped[s].size()});
			tasks.push_back({b, s, std::move(grouped[s])});
			grouped[s].clear();
		}
	}
	PartitionWorkPool& pool = coarseWorkPool_ ? *coarseWorkPool_ : *orderedPool_;
	// A shard's jobs run in enqueue order, i.e. block order; different shards
	// (independent keys) run concurrently on the workers waiting here.
	pool.run(work, [&](size_t index) {
		const Task& task = tasks[index];
		auto& shard = *shards_[task.shard];
		std::lock_guard<std::mutex> guard(shard.mutex);
		for (size_t entry : task.indices) {
			apply_pending_entry(blocks[task.block].entries[entry], shard);
		}
	}, [&submitLock] { submitLock.unlock(); });
	if (submitLock.owns_lock()) { submitLock.unlock(); }
	// This thread may have taken over blocks sealed by other producers. Freeing
	// their reads here would lock those producers' malloc arenas (and did, for up
	// to half of all arena lock waits): hand them back to be freed by their owner.
	for (auto& block : blocks) {
		if (block.owner != self) { return_reads(block.owner, std::move(block.entries)); }
	}
}

void Dereplicate::drain_sealed_blocks() {
	// Caller holds LifecyclePause. Only blocks after a missing (failed) block
	// can still be waiting here; apply them in order.
	std::vector<SealedBlock> ready;
	std::unique_lock<std::mutex> submitLock(submitMutex_, std::defer_lock);
	bool commit = false;
	{
		std::lock_guard<std::mutex> lock(sequenceMutex_);
		if (sealed_.empty()) { return; }
		for (auto& entry : sealed_) { ready.push_back(std::move(entry.second)); }
		nextSubmit_ = sealed_.rbegin()->first + 1;
		sealed_.clear();
		if (parallelGlobalCoarse_) { commit = queue_global_commits(ready); }
		else { submitLock.lock(); }
	}
	if (!parallelGlobalCoarse_) { apply_blocks_ordered(ready, submitLock, nullptr); }
	// No producer is active, so no committer is running either.
	else if (commit) { commit_global_blocks(); }
}

void Dereplicate::BCnamesAdding(Filters* fil) {
	vector<string> refSID = fil->SampleID;
	for (size_t i = 0; i < refSID.size(); i++) {
		barcode_number_to_sample_id_.push_back(refSID[i]);
	}
	//this will in the end be used to synchronize the different barcodes between samples
	curBCoffset = (int)barcode_number_to_sample_id_.size();
	//cerr << curBCoffset << " Added BC list\n";
}
void Dereplicate::reset() {
	LifecyclePause pause(*this);
	flush_all_batches();
	if (merger) { merger->restStats(); }
//	for (size_t i = 0; i < Dnas.size(); i++) { delete Dnas[i]; }
//	Dnas.resize(0);
	//passedSize = 0; 
	tmpCnt = 0;
	//barcode_number_to_sample_id_.resize(0);
	clear_shards();
    if (coarseWorkPool_) { coarseWorkPool_->resetStats(); }
	srch_seq_total_len_.store(0, std::memory_order_release);
	srch_seq_count_.store(0, std::memory_order_release);
    skippedEmptyMates_.store(0, std::memory_order_relaxed);
    skippedEmptyKeys_.store(0, std::memory_order_relaxed);
	lifecycle_pause_.store(false, std::memory_order_release);
	lifecycle_wait_cv_.notify_all();
}
void Dereplicate::clear_shards() {
    // Caller holds LifecyclePause. Freeing the dereplication tables is a large
    // part of the end of a run; shards are independent, so free them
    // concurrently on the otherwise idle workers.
    deferredEvidenceBytes_.store(0, std::memory_order_relaxed);
    auto clearOne = [this](size_t i) {
        auto& shard = *shards_[i];
        std::lock_guard<std::mutex> guard(shard.mutex);
        shard.tracker.clear();
        shard.rejected.clear();
        shard.rejectedTallies.clear();
        shard.deferredCounts.clear();
        shard.deferredEvidence.clear();
        if (shard.coarse) { shard.coarse->clear(); }
        shard.admission = {};
    };
    ThreadPool* pool = nullptr;
    if (shards_.size() > 1 && ThreadPool::configured_thread_count() > 1) {
        ThreadPool& candidate = ThreadPool::instance();
        if (candidate.worker_index_of_current_thread() < 0) { pool = &candidate; }
    }
    if (!pool) {
        for (size_t i = 0; i < shards_.size(); ++i) { clearOne(i); }
        return;
    }
    std::vector<std::future<void>> done;
    done.reserve(shards_.size());
    std::exception_ptr failure;
    for (size_t i = 0; i < shards_.size(); ++i) {
        try { done.push_back(pool->submit([&clearOne, i]() { clearOne(i); })); }
        catch (...) { failure = std::current_exception(); break; }
    }
    for (auto& task : done) {
        try { task.get(); }
        catch (...) { if (!failure) { failure = std::current_exception(); } }
    }
    if (failure) { std::rethrow_exception(failure); }
}

bool DNAuPointerCompare(shared_ptr<DNAunique> l, shared_ptr<DNAunique> r) { 
	return l->totalSum() > r->totalSum(); 
}

std::vector<Dereplicate::DerepOutputGroup> Dereplicate::outputGroups(std::vector<uint64_t>& recoveredCounts, ThreadPool* pool) {
    // Caller holds LifecyclePause and has flushed every producer batch.
    // Every shard fills its own slice: coarse groups of all shards first, then
    // the exact trackers, each in shard order (the order of a serial pass).
    const size_t shardCount = shards_.size();
    std::vector<size_t> coarseStart(shardCount + 1, 0), trackerStart(shardCount + 1, 0);
    size_t total = 0;
    for (size_t s = 0; s < shardCount; ++s) {
        coarseStart[s] = total;
        if (coarse_) {
            const auto& coarse = *shards_[s]->coarse;
            total += coarse.outputsExact() ? coarse.exactGroups().size() : coarse.groups().size();
        }
    }
    coarseStart[shardCount] = total;
    for (size_t s = 0; s < shardCount; ++s) { trackerStart[s] = total; total += shards_[s]->tracker.size(); }
    trackerStart[shardCount] = total;
    vector<DerepOutputGroup> groups(total);
    auto collect = [&](size_t s) {
        const auto& shard = *shards_[s];
        size_t out = coarseStart[s];
        if (coarse_) {
            if (shard.coarse->outputsExact()) {
                for (const auto& [key, assignment] : shard.coarse->exactGroups()) {
                    groups[out++] = {assignment.representatives->totalSum(), key, assignment.representatives.get(), {}};
                }
            } else {
                for (const auto& group : shard.coarse->groups()) {
                    groups[out++] = {group->representatives.totalSum(), group->key(), &group->representatives, {}};
                }
            }
        }
        out = trackerStart[s];
        for (const auto& entry : shard.tracker) {
            assert(entry.second);
            groups[out++] = {entry.second->totalSum(), entry.first, entry.second.get(), {}};
        }
    };
    if (pool && shardCount > 1) { runOnPool(*pool, shardCount, collect); }
    else { for (size_t s = 0; s < shardCount; ++s) { collect(s); } }
    auto addCount = [](uint64_t& target, uint64_t count) {
        if (count > std::numeric_limits<uint64_t>::max() - target) {
            throw std::overflow_error("Dereplication recovered count overflow");
        }
        target += count;
    };
    // Deferred non-passing reads join the group they are attached to: counts
    // now, quality/variant evidence (members) and sample tallies at write time.
    auto recoverMembers = [&addCount](DerepOutputGroup& group, DNAuniqSet* deferred, std::vector<uint64_t>& recovered) {
        addCount(group.count, deferred->totalSum());
        group.members.push_back(deferred);
        for (auto subgroup = deferred->begin(); subgroup != deferred->end(); ++subgroup) {
            for (const auto& occurrence : subgroup->second->getDerepMap()) {
                if (occurrence.first >= 0 && static_cast<size_t>(occurrence.first) < recovered.size()) {
                    addCount(recovered[occurrence.first], occurrence.second);
                }
            }
        }
    };
    auto recoverCounts = [&addCount](DerepOutputGroup& group, int sample, uint64_t copies, std::vector<uint64_t>& recovered) {
        addCount(group.count, copies);
        group.recovered.emplace_back(sample, copies);
        if (sample >= 0 && static_cast<size_t>(sample) < recovered.size()) { addCount(recovered[sample], copies); }
    };
    if (!prefixMode_) {
        if (!defersNonPassing()) { return groups; }
        // Exact keys: a deferred read belongs to the passing group of its own
        // key, which lives in the same shard (same key hash). Keys without a
        // passing group are dropped here. Shards touch disjoint group ranges.
        std::vector<std::vector<uint64_t>> recoveredPerShard(shardCount, std::vector<uint64_t>(recoveredCounts.size(), 0));
        auto attach = [&](size_t s) {
            const auto& deferredCounts = shards_[s]->deferredCounts;
            const auto& deferredEvidence = shards_[s]->deferredEvidence;
            if (deferredCounts.size() == 0 && deferredEvidence.empty()) { return; }
            auto recoverGroup = [&](DerepOutputGroup& group) {
                const uint64_t key = deferredKeyHash(group.key);
                if (auto reads = deferredEvidence.find(key); reads != deferredEvidence.end()) {
                    for (const auto& read : reads->second) {
                        addCount(group.count, 1);
                        group.evidence.push_back(&read);
                        if (read.sample >= 0 && static_cast<size_t>(read.sample) < recoveredPerShard[s].size()) {
                            addCount(recoveredPerShard[s][read.sample], 1);
                        }
                    }
                }
                if (auto one = deferredCounts.single.find(key); one != deferredCounts.single.end()) {
                    recoverCounts(group, one->second.sample, one->second.copies, recoveredPerShard[s]);
                }
                if (auto several = deferredCounts.multi.find(key); several != deferredCounts.multi.end()) {
                    for (const auto& [sample, copies] : several->second) { recoverCounts(group, sample, copies, recoveredPerShard[s]); }
                }
            };
            for (size_t g = coarseStart[s]; g < coarseStart[s + 1]; ++g) { recoverGroup(groups[g]); }
            for (size_t g = trackerStart[s]; g < trackerStart[s + 1]; ++g) { recoverGroup(groups[g]); }
        };
        if (pool && shardCount > 1) { runOnPool(*pool, shardCount, attach); }
        else { for (size_t s = 0; s < shardCount; ++s) { attach(s); } }
        for (const auto& part : recoveredPerShard) {
            for (size_t j = 0; j < recoveredCounts.size(); ++j) { addCount(recoveredCounts[j], part[j]); }
        }
        return groups;
    }
    auto keyLess = [](const auto& a, const auto& b) { return a.key < b.key; };
    sortOnPool(groups, keyLess, pool); // keys are unique: a total order
    vector<DerepOutputGroup> anchors;
    anchors.reserve(groups.size());
    // A prefix is adjacent to its first extension. Traverse backwards so its
    // assignment already names a maximal passing sequence. Conflicting longer
    // sequences stay separate; an ambiguous prefix goes to the first in
    // lexicographic order, independently of read arrival or worker scheduling.
    for (size_t i = groups.size(); i-- > 0;) {
        if (i + 1 < groups.size() && groups[i + 1].key.starts_with(groups[i].key)) {
            auto& anchor = anchors.back();
            addCount(anchor.count, groups[i].count);
            anchor.members.push_back(groups[i].representative);
        } else {
            anchors.push_back(std::move(groups[i]));
        }
    }
    // Anchors were collected walking the key-sorted groups backwards, so they
    // are in descending key order; with unique keys, reversing equals sorting.
    std::reverse(anchors.begin(), anchors.end());
    assert(std::is_sorted(anchors.begin(), anchors.end(), keyLess));
    // A deferred key attaches to the anchor it is a prefix of (the first in
    // lexicographic order), or to the anchor that is its prefix; else it is dropped.
    auto findAnchor = [&anchors](std::string_view key) {
        auto match = lower_bound(anchors.begin(), anchors.end(), key,
            [](const auto& group, std::string_view value) { return group.key < value; });
        if (match == anchors.end() || !match->key.starts_with(key)) {
            if (match == anchors.begin()) { return anchors.end(); }
            --match;
            if (!key.starts_with(match->key)) { return anchors.end(); }
        }
        return match;
    };
    for (const auto& shard : shards_) {
        for (const auto& entry : shard->rejected) {
            auto match = findAnchor(entry.first);
            if (match != anchors.end()) { recoverMembers(*match, entry.second.get(), recoveredCounts); }
        }
        for (const auto& entry : shard->rejectedTallies) {
            auto match = findAnchor(entry.first);
            if (match == anchors.end()) { continue; }
            for (const auto& [sample, copies] : entry.second.counts) { recoverCounts(*match, sample, copies, recoveredCounts); }
        }
    }
    return anchors;
}

void Dereplicate::finishMap() {
	LifecyclePause pause(*this);
	flush_all_batches();
	//at this point we can onl be sure that barcode_number_to_sample_id_ is finished
	//hence now it the point to add this to map
	//passedSize = 0; 
	tmpCnt = 0;
	clear_shards();
    if (coarseWorkPool_) { coarseWorkPool_->resetStats(); }
	srch_seq_total_len_.store(0, std::memory_order_release);
	srch_seq_count_.store(0, std::memory_order_release);
    skippedEmptyMates_.store(0, std::memory_order_relaxed);
    skippedEmptyKeys_.store(0, std::memory_order_relaxed);
	const string tmpMap = mapF + ".tmp";
	const string bakMap = mapF + ".bak";
	std::ifstream inputFile(mapF.c_str(), ios::in);
	if (!inputFile.is_open()) { throw std::runtime_error("Could not open dereplication map: " + mapF); }
    inputFile.exceptions(ios::badbit);
	std::ofstream outputFile(tmpMap.c_str(), ios::out | ios::trunc);
    if (!outputFile.is_open()) { throw std::runtime_error("Could not create temporary dereplication map: " + tmpMap); }

	outputFile << "#SMPLS";
	bool bCombiSmpl = mainFilter->combineSamples();
	vector<int> smplId2comb(0, 0);
	unordered_map<string, int> & combiMapCollectGrp = mainFilter->combiMapCollectGrp;
	if (!bCombiSmpl) {
    for (size_t i = 0; i < barcode_number_to_sample_id_.size(); i++) {
			outputFile << "\t" + std::to_string((int)i) + ":" + barcode_number_to_sample_id_[i];
		}
	}
	else {
		smplId2comb = mainFilter->combiSmplConvergeVec(barcode_number_to_sample_id_);
		if (barcode_number_to_sample_id_.size() != smplId2comb.size()) {
			cerr << "FATAL: barcode_number_to_sample_id_ != smplId2comb\n"; exit(234);
		}
        for (auto IT = combiMapCollectGrp.begin(); IT != combiMapCollectGrp.end(); IT++) {
			outputFile << "\t" << std::to_string(IT->second) + ":" + IT->first;
		}
	}
	outputFile << "\n";
    // Copy explicitly so an empty map is valid and read errors are checked.
    char buffer[64 * 1024];
    while (inputFile.read(buffer, sizeof(buffer)) || inputFile.gcount() > 0) {
        outputFile.write(buffer, inputFile.gcount());
    }
    if (!inputFile.eof()) { throw std::runtime_error("Could not read dereplication map: " + mapF); }
    inputFile.close();
    finalizeOutputFile(outputFile, tmpMap);

    // Never promote an incomplete file or discard the original on rename failure.
    if (std::rename(mapF.c_str(), bakMap.c_str()) != 0) {
        throw std::runtime_error("Could not back up dereplication map: " + mapF);
    }
    if (std::rename(tmpMap.c_str(), mapF.c_str()) != 0) {
        const bool restored = std::rename(bakMap.c_str(), mapF.c_str()) == 0;
        throw std::runtime_error("Could not promote temporary dereplication map: " + tmpMap +
            (restored ? "; original map restored" : "; original map retained at " + bakMap));
    }
    std::remove(bakMap.c_str());
    b_wroteMapHD = true;

	return;
}
string Dereplicate::writeDereplDNA(Filters* mf, string SRblock) {
	if (outfile.empty() || mapF.empty()) {
		throw std::runtime_error("Dereplicate output paths are not initialized before writeDereplDNA().");
	}

	LifecyclePause pause(*this);
	flush_all_batches();

	ofstream of, omaps, of2, ofRest, of2p2, of_merged;
	cerr << "\nEvaluating and writing dereplicated reads..\n";
	int fastqVer = mf->getuserReqFastqOutVer();
	string srBlockSafe = SRblock;
	trim(srBlockSafe);
	if (srBlockSafe.size() >= 3 &&
		static_cast<unsigned char>(srBlockSafe[0]) == 0xEF &&
		static_cast<unsigned char>(srBlockSafe[1]) == 0xBB &&
		static_cast<unsigned char>(srBlockSafe[2]) == 0xBF) {
		srBlockSafe.erase(0, 3);
	}
	srBlockSafe.erase(std::remove_if(srBlockSafe.begin(), srBlockSafe.end(), [](char c) {
		unsigned char uc = static_cast<unsigned char>(c);
		if (std::isspace(uc)) {
			return true;
		}
		switch (c) {
		case '/': case '\\': case ':': case '*': case '?': case '"': case '<': case '>': case '|':
			return true;
		default:
			return false;
		}
	}), srBlockSafe.end());
	
	//set the correct file names for dereplication output files
	string baseOF = std::filesystem::path(outfile).replace_extension().string();
	string baseOF2 = baseOF;	
	if (b_derepPerSR && srBlockSafe != "") {//if writing out per SRblock, needs different filename for map and 
		baseOF2 += "." + srBlockSafe;
	}	
	string fileSuff = std::filesystem::path(outfile).extension().string();
	string outfile2 = baseOF2 + fileSuff;
    string outfile_merged = baseOF2 + ".merg" + fileSuff;
	auto ensureOpen = [](const std::ofstream& stream, const std::string& path, const char* label) {
		if (!stream.is_open() || stream.fail()) {
			throw std::runtime_error(std::string("Failed to open ") + label + " file: " + path);
		}
	};

    std::ofstream diffOutput, subclusterOutput, subclusterOutputMate;
    const bool variantHq = coarse_ && coarse_->storesQualities();
    const bool pairedVariants = variantHq && b_pairedInput;
    const string subclusterPath = baseOF2 + (pairedVariants ? ".subclusters.1.fna" : ".subclusters.fna");
    const string subclusterMatePath = baseOF2 + ".subclusters.2.fna";
    const string diffPath = baseOF2 + ".diff";
    if (coarse_ && coarse_->storesDiffs()) {
        diffOutput.open(diffPath, ios::out | ios::binary | ios::trunc);
        ensureOpen(diffOutput, diffPath, "nucleotide differences");
        if (writeSubclusters_) {
            subclusterOutput.open(subclusterPath, ios::out | ios::trunc);
            ensureOpen(subclusterOutput, subclusterPath, "exact subcluster FASTA");
            if (pairedVariants) {
                subclusterOutputMate.open(subclusterMatePath, ios::out | ios::trunc);
                ensureOpen(subclusterOutputMate, subclusterMatePath, "exact subcluster mate FASTA");
            }
        }
    }
	// OPEN OUTPUT FILE STREAMS
	const auto cumulativeMode = b_derepPerSR ? ios::app : (ios::out | ios::trunc);
	omaps.open(mapF.c_str(), cumulativeMode);//| ios::binary
	ensureOpen(omaps, mapF, "derep map");
	//actual dereplicated reads as required by clustering algo (uparse, cdhit, dada2,...)
	of.open(outfile2.c_str(), ios::out);
	ensureOpen(of, outfile2, "dereplicated output");
	//reads that did not pass dereplication min number
	ofRest.open(outRest.c_str(), cumulativeMode);
	ensureOpen(ofRest, outRest, "dereplicated rest");
	//copy of high quality read pairs -> needed for Seed extension step
	of2.open(outHQf.c_str(), cumulativeMode);
	ensureOpen(of2, outHQf, "high-quality output");
	const bool can_merge_derep = b_merge_pairs_derep_ && (merger != nullptr);
	if (b_merge_pairs_derep_ && merger == nullptr) {
		cerr << "Warning: -merge_pairs_derep enabled but merger not initialized; continuing without merge for derep output.\n";
	}
	if (can_merge_derep) {
        of_merged.open(outfile_merged.c_str(), ios::out);
		ensureOpen(of_merged, outfile_merged, "merged dereplicated output");
    }
	if (b_pairedInput) {
		of2p2.open(outHQf_p2, cumulativeMode);
		ensureOpen(of2p2, outHQf_p2, "high-quality output read2");
	}



	//sample specific derep filter
	bool dereplicate_sample_specific(false);
	vector<int> dereplicate_sample_specific_indices;
	{
		// derepMinNum follows mapping rows, but sample counts use dereplicator ids
		// (file processing order): translate by sample name.
		const vector<int> byMapRow = mf->getDrerepSampleSpecifity();
		if (!byMapRow.empty()) {
			std::unordered_map<string, int> rowOfSample;
			for (size_t i = 0; i < mf->SampleID.size() && i < byMapRow.size(); ++i) {
				rowOfSample.emplace(mf->SampleID[i], byMapRow[i]);
			}
			dereplicate_sample_specific_indices.assign(barcode_number_to_sample_id_.size(), -1);
			for (size_t j = 0; j < barcode_number_to_sample_id_.size(); ++j) {
				auto row = rowOfSample.find(barcode_number_to_sample_id_[j]);
				if (row != rowOfSample.end()) { dereplicate_sample_specific_indices[j] = row->second; }
			}
		}
	}
	if (!dereplicate_sample_specific_indices.empty()) { dereplicate_sample_specific = true; }

	//combiner relevant vars
	bool bCombiSmpl = mf->combineSamples();
	vector<int> smplId2comb(0,0);
	unordered_map<string, int> & combiMapCollectGrp = mf->combiMapCollectGrp;
	if (bCombiSmpl) {
		smplId2comb = mf->combiSmplConvergeVec(barcode_number_to_sample_id_);
	}
	//vector<string> refSID = mf->SampleID;
	if (false && !b_wroteMapHD) {
		omaps << "#SMPLS";
		if (!bCombiSmpl) {
			for (size_t i = 0; i < barcode_number_to_sample_id_.size(); i++) {
				omaps << "\t" + itos((int)i) + ":" + barcode_number_to_sample_id_[i];
			}
		}
		else {
			if (barcode_number_to_sample_id_.size() != smplId2comb.size()) {
				cerr << "FATAL: barcode_number_to_sample_id_ != smplId2comb\n"; exit(234);
			}
        for (auto IT = combiMapCollectGrp.begin(); IT != combiMapCollectGrp.end(); IT++) {
			omaps << "\t" << std::to_string(IT->second) + ":" + IT->first;
			}
		}
		omaps << "\n";
		b_wroteMapHD = true;
	}
	
	//convert to vector, that can than be written out
	// Sort small references; materialize only one independent output group at a time.
    // Collecting and sorting the output groups may use the pool (never from one of
    // its own workers).
    ThreadPool* sortPool = nullptr;
    if (ThreadPool::configured_thread_count() > 1) {
        ThreadPool& pool = ThreadPool::instance();
        if (pool.worker_index_of_current_thread() < 0) { sortPool = &pool; }
    }
    vector<uint64_t> recoveredCounts(barcode_number_to_sample_id_.size(), 0);
    CoarseDereplication::Refinement refinement;
    std::vector<DerepOutputGroup> dereplicated_dnas;
    if (coarse_ && coarse_->reassigns()) {
        refinement = coarse_->refine(searchWithMerg, searchLength);
        for (size_t i = 0; i < refinement.representatives.size(); ++i) {
            auto& dna = refinement.representatives[i];
            if (dna) { dereplicated_dnas.push_back({dna->totalSum(), coarse_->groups()[i]->key(), nullptr, {}, dna}); }
        }
    } else { dereplicated_dnas = outputGroups(recoveredCounts, sortPool); }
    // Keys are unique, so this order is total and independent of how the sort is split.
    sortOnPool(dereplicated_dnas, [](const auto& left, const auto& right) {
        return left.count != right.count ? left.count > right.count : left.key < right.key;
    }, sortPool);
	size_t passedSize = 0; size_t notPassedSize = 0; size_t passed_hits(0);
	size_t passedMerge(0);
	//bool thrHit = false;

	//sanity check
	vector<uint64_t> counts_per_sample(barcode_number_to_sample_id_.size(), 0);
	vector<uint64_t> map_combined_counts;
	vector<int> map_touched_groups;
	if (!smplId2comb.empty()) {
		int max_combined_id = -1;
		for (size_t i = 0; i < smplId2comb.size(); ++i) {
			if (smplId2comb[i] > max_combined_id) {
				max_combined_id = smplId2comb[i];
			}
		}
		if (max_combined_id >= 0) {
			map_combined_counts.assign(static_cast<size_t>(max_combined_id) + 1, 0);
			map_touched_groups.reserve(std::min(static_cast<size_t>(max_combined_id) + 1, smplId2comb.size()));
		}
	}
	int total_count(0);
	int passed_count(0);
	// Groups are independent: format them in chunks, in parallel when workers
	// are available, and write the chunks in order so every output file is
	// byte-identical to formatting them one after another.
	struct ChunkSinks {
		std::ostream *diff, *subcluster, *subclusterMate, *hq1, *hq2, *main, *rest, *merged, *map;
	};
	struct ChunkCounts {
		int total = 0;
		size_t passedHits = 0, passedMerge = 0, passedSize = 0, notPassedSize = 0;
		vector<uint64_t> countsPerSample;
	};
	SegmentPool segments;
	struct ChunkOutput {
		explicit ChunkOutput(SegmentPool& pool) : diff(pool), subcluster(pool), subclusterMate(pool), hq1(pool), hq2(pool),
			main(pool), rest(pool), merged(pool), map(pool) {}
		SegmentStream diff, subcluster, subclusterMate, hq1, hq2, main, rest, merged, map;
		ChunkCounts counts;
		ChunkSinks sinks() {
			return {&diff.stream, &subcluster.stream, &subclusterMate.stream, &hq1.stream, &hq2.stream,
				&main.stream, &rest.stream, &merged.stream, &map.stream};
		}
	};
	std::mutex mergeCheckMutex; // the main filter's merged-read checks keep statistics
	// A large group's exact-variant FASTQ, left for the caller to write in parallel ranges.
	struct DeferredVariants {
		std::shared_ptr<const diffDNA> variants;
		std::string id;
	};
	auto formatChunk = [&](size_t begin, size_t end, const ChunkSinks& out, ChunkCounts& counts, DeferredVariants* defer = nullptr) {
		counts.countsPerSample.assign(counts_per_sample.size(), 0);
		vector<uint64_t> mapCombinedCounts(map_combined_counts.size(), 0);
		vector<int> mapTouchedGroups;
		mapTouchedGroups.reserve(map_touched_groups.capacity());
	for (size_t i = begin; i < end; i++) {
        const auto& group = dereplicated_dnas[i];
        shared_ptr<DNAunique> dna = group.refined ? group.refined : group.representative->best(true, !b_derep_as_fasta_);
        if (!dna) { continue; }
        for (const auto* member : group.members) { member->addEvidenceTo(dna, !b_derep_as_fasta_); }
        for (const auto& [sample, copies] : group.recovered) { dna->addSampleCount(sample, copies); }
        // Deferred reads of variant-retaining modes become exact variants of this
        // parent, exactly as a failed repeat of an established key is recorded.
        for (const auto* read : group.evidence) {
            if (coarse_->storesQualities() && read->dna2) { dna->recordDiffPairQuality(*read->dna, *read->dna2, read->sample); }
            else if (coarse_->storesQualities()) { dna->recordDiffQuality(*read->dna, read->sample); }
            else {
                const auto& source = searchWithMerg && read->merged ? read->merged : read->dna;
                dna->recordDiffDNA(source->getSeqPseudoView(), read->sample);
            }
        }
        counts.total ++;
		// R1's header before the size annotation; stored R2 headers are rebuilt from it.
		const string r1Header = dna->storedHeader();
		dna->Count2Head(b_usearch_fmt);
        if (auto differences = dna->getDiffDNA(); differences && coarse_) {
            if (coarse_->storesDiffs()) {
                differences->write(*out.diff, dna->getId());
                if (writeSubclusters_) {
                    if (pairedVariants) { differences->writeFasta(*out.subcluster, *out.subclusterMate, dna->getId()); }
                    else { differences->writeFasta(*out.subcluster, dna->getId()); }
                }
            }
            if (variantHq && defer) {
                defer->variants = differences;
                defer->id = dna->getId();
            } else if (variantHq) {
                // Use the normal cumulative HQ streams for every retained exact
                // variant. Parent IDs/counts still resolve through the usual map.
                if (pairedVariants) { differences->writeFastq(*out.hq1, *out.hq2, dna->getId(), fastqVer); }
                else { differences->writeFastq(*out.hq1, dna->getId(), fastqVer); }
            }
        }
		shared_ptr<DerepRead> dna_merged = nullptr;
		std::ostream* derepNowOut;
		if ((dereplicate_sample_specific && 
			dna->pass_deprep_smplSpc(dereplicate_sample_specific_indices)) ||
            pass_deprep_conditions(dna) ) {
			//we do have a real derep that needs to be clustered..
			
			counts.passedHits++;
			counts.passedSize += dna->totalSum();


			//only do merge here, because these are known good dereps already
            if (can_merge_derep && dna->merge_seed_pos_ >= 0) {
                dna_merged = dna->getMerge();
                if (!dna_merged) {
                    auto consensus = merger->merge(dna, dna->getPair());
                    if (consensus) { dna_merged = make_shared<DerepRead>(consensus); }
                    if (dna_merged) {
                        std::lock_guard<std::mutex> lock(mergeCheckMutex);
                        if (!mf->checkMergedDNA(dna_merged)) { dna_merged.reset(); }
                    }
                }
                if (dna_merged) {
                    std::lock_guard<std::mutex> lock(mergeCheckMutex);
                    if (!mf->checkMergedLength(dna_merged)) { dna_merged.reset(); }
                }
                if (dna_merged) { dna_merged->setHeader(dna->getId()); dna->adoptMergeStats(*dna_merged); }
			}
			derepNowOut = out.main;
			if (dna->getMergeLength() > 0) {
				counts.passedMerge++;
			}
		} else {
			counts.notPassedSize += dna->totalSum();
			//dna->writeSeq(ofRest, b_singleLine);
			derepNowOut = out.rest;
		}

		if (b_derep_as_fasta_)
			if (dna_merged) {//either or mechanic
				dna_merged->writeSeq_derep(*out.merged, b_singleLine);
			} else {
				dna->writeSeq_derep((*derepNowOut), b_singleLine);
			}
		else {
			if (dna_merged) {
				dna_merged->prepareDerepQualities(fastqVer);
				dna_merged->writeDerepFastQ(*out.merged);
			} else {
				dna->prepareDerepQualities(fastqVer);
				dna->writeDerepFastQ((*derepNowOut));
			}

		}


		//map is written no matter what
		dna->writeMap(*out.map, dna->getId(), counts.countsPerSample, smplId2comb, mapCombinedCounts, mapTouchedGroups);

        if (!variantHq) {
            // HQ mates are full-length intermediates for seed extension. Preserve
            // their overlap even when clustering uses truncated reads; physical
            // barcode/primer cuts remain in the backing sequence.
            dna->resetTruncation();
            dna->prepareWrite(fastqVer);
            dna->writeFastQ(*out.hq1);
            dna->releaseRawQualitiesIfPromoted();

            if (b_pairedInput) {
                shared_ptr<DerepRead> oD = dna->getPair();
                if (oD != nullptr) {
                    oD->restoreMateHeader(r1Header);
                    oD->resetTruncation();
                    oD->prepareWrite(fastqVer);
                    oD->writeFastQ(*out.hq2);
                    oD->releaseRawQualitiesIfPromoted();
                } else {
                    shared_ptr<DNA> tmp = make_shared< DNA>("", dna->getId());
                    tmp->writeFastQEmpty(*out.hq2);

                }
            }
        }
		
	}
	};
	auto addCounts = [&](const ChunkCounts& counts) {
		total_count += counts.total;
		passed_hits += counts.passedHits;
		passedMerge += counts.passedMerge;
		passedSize += counts.passedSize;
		notPassedSize += counts.notPassedSize;
		for (size_t j = 0; j < counts_per_sample.size(); ++j) { counts_per_sample[j] += counts.countsPerSample[j]; }
	};
	auto writeChunk = [&](ChunkOutput& out) {
		auto put = [](std::ostream& target, SegmentStream& source) { source.buffer.writeTo(target); };
		if (diffOutput.is_open()) { put(diffOutput, out.diff); }
		if (subclusterOutput.is_open()) { put(subclusterOutput, out.subcluster); }
		if (subclusterOutputMate.is_open()) { put(subclusterOutputMate, out.subclusterMate); }
		put(of2, out.hq1);
		if (b_pairedInput) { put(of2p2, out.hq2); }
		put(of, out.main);
		put(ofRest, out.rest);
		if (of_merged.is_open()) { put(of_merged, out.merged); }
		put(omaps, out.map);
		addCounts(out.counts);
	};
	// Record writers consult the runtime options; refresh their lazily built
	// cache here so that concurrent formatting only reads it.
	if (OptContainer* options = getGlobalOptContainer()) { (void)options->runtime(); }
	ThreadPool* formatPool = nullptr;
	// Variant export writes (and, to sort them, decodes) every retained exact
	// variant of a group; see the chunk plan below for how that stays bounded.
	const bool exportsVariants = coarse_ && (coarse_->storesDiffs() || coarse_->storesQualities());
	if (dereplicated_dnas.size() > 1 && ThreadPool::configured_thread_count() > 1) {
		ThreadPool& pool = ThreadPool::instance();
		// Never wait for the pool from inside one of its own workers.
		if (pool.worker_index_of_current_thread() < 0) { formatPool = &pool; }
	}
	if (!formatPool) {
		// One worker: write straight into the files, without buffering.
		ChunkSinks direct{&diffOutput, &subclusterOutput, &subclusterOutputMate, &of2, &of2p2, &of, &ofRest, &of_merged, &omaps};
		ChunkCounts counts;
		formatChunk(0, dereplicated_dnas.size(), direct, counts);
		addCounts(counts);
	} else {
		// Chunks are formatted on the pool and written in order; a small window
		// of chunks keeps buffered output bounded.
		const size_t window = formatPool->size() + 2;
		struct Chunk { size_t begin, end; bool large; };
		std::vector<Chunk> chunks;
		constexpr size_t kMaxGroupsPerChunk = 512;
		if (exportsVariants) {
			// A group's variant output, and the variants decoded to sort it,
			// grow with its retained variants, which never exceed its reads.
			// Chunks hold at most kChunkReads reads, so neither the buffered
			// chunks nor any worker's formatting memory grow with the dataset.
			// Larger groups are written by this thread one at a time, as in
			// the serial path; their sorted variants are formatted in
			// kChunkReads ranges on the pool.
			constexpr uint64_t kChunkReads = 1024;
			size_t begin = 0;
			uint64_t reads = 0;
			for (size_t i = 0; i < dereplicated_dnas.size(); ++i) {
				const uint64_t count = dereplicated_dnas[i].count;
				if (begin < i && (count > kChunkReads || reads + count > kChunkReads || i - begin >= kMaxGroupsPerChunk)) {
					chunks.push_back({begin, i, false});
					begin = i;
					reads = 0;
				}
				if (count > kChunkReads) {
					chunks.push_back({i, i + 1, true});
					begin = i + 1;
					continue;
				}
				reads += count;
			}
			if (begin < dereplicated_dnas.size()) { chunks.push_back({begin, dereplicated_dnas.size(), false}); }
		} else {
			// Chunks of similar abundance: output size follows the number of
			// retained reads, and groups are sorted by abundance, so fixed group
			// counts would put most of the output into the first chunk.
			uint64_t totalAbundance = 0;
			for (const auto& group : dereplicated_dnas) { totalAbundance += group.count; }
			const uint64_t chunkAbundance = std::max<uint64_t>(1, totalAbundance / (32 * formatPool->size()));
			size_t begin = 0;
			uint64_t accumulated = 0;
			for (size_t i = 0; i < dereplicated_dnas.size(); ++i) {
				accumulated += dereplicated_dnas[i].count;
				if (accumulated >= chunkAbundance || i + 1 - begin >= kMaxGroupsPerChunk) {
					chunks.push_back({begin, i + 1, false});
					begin = i + 1;
					accumulated = 0;
				}
			}
			if (begin < dereplicated_dnas.size()) { chunks.push_back({begin, dereplicated_dnas.size(), false}); }
		}
		// Waits for every task in `pending` (they reference this frame), then rethrows `failure`.
		auto settle = [](auto& pending, std::exception_ptr failure) {
			for (auto& entry : pending) {
				try { if (entry.second.valid()) { entry.second.get(); } }
				catch (...) {}
			}
			if (failure) { std::rethrow_exception(failure); }
		};
		// Decoding and sorting a large group's variants run as pool tasks too.
		const diffDNA::ParallelFor parallelFor = [&](size_t count, const std::function<void(size_t)>& task) {
			std::vector<std::pair<std::nullptr_t, std::future<void>>> done;
			done.reserve(count);
			std::exception_ptr failure;
			try {
				for (size_t i = 0; i < count; ++i) { done.emplace_back(nullptr, formatPool->submit([&task, i]() { task(i); })); }
				for (auto& entry : done) { entry.second.get(); }
			}
			catch (...) { failure = std::current_exception(); }
			settle(done, failure);
		};
		// The one large group being written: its sorted variants are shared by
		// the range tasks and released before the next large group is decoded.
		auto writeVariantRanges = [&](const diffDNA& variants, const std::string& id) {
			const auto fragments = variants.sortedFragments(&parallelFor);
			constexpr size_t kRange = 1024;
			struct RangeOutput {
				explicit RangeOutput(SegmentPool& pool) : hq1(pool), hq2(pool) {}
				SegmentStream hq1, hq2;
			};
			std::deque<std::pair<std::unique_ptr<RangeOutput>, std::future<void>>> pending;
			size_t next = 0;
			auto submitRange = [&]() {
				auto out = std::make_unique<RangeOutput>(segments);
				RangeOutput* target = out.get();
				const size_t begin = next, end = std::min(fragments.size(), next + kRange);
				auto done = formatPool->submit([&variants, &fragments, &id, target, begin, end, pairedVariants, fastqVer]() {
					if (pairedVariants) { variants.writeFastq(target->hq1.stream, target->hq2.stream, id, fragments, begin, end, fastqVer); }
					else { variants.writeFastq(target->hq1.stream, id, fragments, begin, end, fastqVer); }
				});
				pending.emplace_back(std::move(out), std::move(done));
				next = end;
			};
			std::exception_ptr failure;
			try {
				while (next < fragments.size() && pending.size() < window) { submitRange(); }
				while (!pending.empty()) {
					pending.front().second.get();
					pending.front().first->hq1.buffer.writeTo(of2);
					if (pairedVariants) { pending.front().first->hq2.buffer.writeTo(of2p2); }
					pending.pop_front();
					if (next < fragments.size()) { submitRange(); }
				}
			}
			catch (...) { failure = std::current_exception(); }
			settle(pending, failure);
		};
		auto writeLargeGroup = [&](size_t index) {
			// Every earlier chunk is written: write straight into the files.
			ChunkSinks direct{&diffOutput, &subclusterOutput, &subclusterOutputMate, &of2, &of2p2, &of, &ofRest, &of_merged, &omaps};
			ChunkCounts counts;
			DeferredVariants deferred;
			formatChunk(index, index + 1, direct, counts, &deferred);
			addCounts(counts);
			if (deferred.variants) { writeVariantRanges(*deferred.variants, deferred.id); }
		};
		// Large chunks hold a window place but are written by this thread.
		std::deque<std::pair<std::unique_ptr<ChunkOutput>, std::future<void>>> inflight;
		size_t nextChunk = 0;
		auto submitNext = [&]() {
			const Chunk chunk = chunks[nextChunk++];
			if (chunk.large) { inflight.emplace_back(nullptr, std::future<void>{}); return; }
			auto out = std::make_unique<ChunkOutput>(segments);
			ChunkOutput* target = out.get();
			auto done = formatPool->submit([&formatChunk, target, chunk]() {
				formatChunk(chunk.begin, chunk.end, target->sinks(), target->counts);
			});
			inflight.emplace_back(std::move(out), std::move(done));
		};
		std::exception_ptr failure;
		try {
			size_t written = 0;
			while (nextChunk < chunks.size() && inflight.size() < window) { submitNext(); }
			while (!inflight.empty()) {
				if (inflight.front().first) {
					inflight.front().second.get();
					writeChunk(*inflight.front().first);
				} else {
					writeLargeGroup(chunks[written].begin);
				}
				inflight.pop_front();
				++written;
				if (nextChunk < chunks.size()) { submitNext(); }
			}
		}
		catch (...) {
			failure = std::current_exception();
		}
		settle(inflight, failure);
	}
	if (of.fail() || omaps.fail() || ofRest.fail() || of2.fail() || (b_pairedInput && of2p2.fail()) || (can_merge_derep && of_merged.fail())) {
		throw std::runtime_error("I/O failure while writing dereplicated outputs.");
	}
//	if (tmpCnt != passedSize) {
//		cerr << "Counting failed\n" << tmpCnt << " " << passedSize<<endl;
//	}
	if (diffOutput.is_open()) { finalizeOutputFile(diffOutput, diffPath); }
    if (subclusterOutput.is_open()) { finalizeOutputFile(subclusterOutput, subclusterPath); }
    if (subclusterOutputMate.is_open()) { finalizeOutputFile(subclusterOutputMate, subclusterMatePath); }
	finalizeOutputFile(of, outfile2);
    finalizeOutputFile(omaps, mapF);
    finalizeOutputFile(ofRest, outRest);
    finalizeOutputFile(of2, outHQf);
    if (b_pairedInput) { finalizeOutputFile(of2p2, outHQf_p2); }
    if (can_merge_derep) { finalizeOutputFile(of_merged, outfile_merged); }
	float avgSize = passed_hits > 0 ? (float)passedSize / (float)(passed_hits) : 0.f;
    string report = "";
	auto formatCount64 = [](uint64_t value) {
		std::string s = std::to_string(value);
		for (int pos = static_cast<int>(s.size()) - 3; pos > 0; pos -= 3) {
			s.insert(static_cast<size_t>(pos), ",");
		}
		return s;
	};
	string N_notPassed = intwithcommas(int(dereplicated_dnas.size() - passed_hits));
	string N_total = intwithcommas(int(dereplicated_dnas.size()));
	string N_passed = intwithcommas((int)passed_hits);
	uint64_t srchSeqCount = srch_seq_count_.load(std::memory_order_acquire);
	uint64_t srchSeqTotalLen = srch_seq_total_len_.load(std::memory_order_acquire);
	double avgSrchSeqLen = srchSeqCount > 0 ? static_cast<double>(srchSeqTotalLen) / static_cast<double>(srchSeqCount) : 0.0;
	std::ostringstream avgSrchSeqLenStream;
	avgSrchSeqLenStream.setf(std::ios::fixed);
	avgSrchSeqLenStream.precision(2);
	avgSrchSeqLenStream << avgSrchSeqLen;
	report += "Dereplication: " + N_passed + " unique sequences (avg size ";
	report += intwithcommas((int)avgSize) + "; " + formatCount64(static_cast<uint64_t>(passedSize)) + " counts, ";
	report += formatCount64(static_cast<uint64_t>(passedMerge)) + " merged)\n";
    report += "Dereplication abundance: " + formatCount64(passedSize + notPassedSize) +
        " total in map; " + formatCount64(passedSize) + " passing; " + formatCount64(notPassedSize) + " in rest\n";
    if (coarse_ && coarse_->storesQualities()) {
        report += "Retention exclusions: " + formatCount64(skippedEmptyMates_.load(std::memory_order_relaxed)) +
            " pairs with empty R2 skipped before admission\n";
    }
    Shard::AdmissionStats admission;
    for (const auto& shard : shards_) {
        const auto& local = shard->admission;
        admission.high += local.high; admission.mid += local.mid; admission.failed += local.failed;
        admission.unmatched += local.unmatched; admission.deferred += local.deferred; admission.mergedKeys += local.mergedKeys;
        admission.published += local.published; admission.deferredDropped += local.deferredDropped;
    }
    report += "Dereplication search policy: " + string(searchWithMerg ? "merged when available" : "R1 only") +
        "; " + (coarse_ && !coarse_->outputsExact() ? "coarse clusters" : prefixMode_ ? "exact prefix consolidation" : "equal exact keys") +
        "; search cap " + to_string(searchLength) + "\n";
    report += "Dereplication key sources: " + formatCount64(admission.mergedKeys) + " merged; " +
        formatCount64(srchSeqCount - admission.mergedKeys) + " R1; " +
        formatCount64(skippedEmptyKeys_.load(std::memory_order_relaxed)) + " empty keys skipped\n";
    report += "Dereplication admission before prefix recovery: " + formatCount64(admission.high) + " high-quality R1; " +
        formatCount64(admission.mid) + " mid-quality R1; " + formatCount64(admission.failed) + " failed R1; " +
        formatCount64(admission.unmatched) + " unmatched; " + formatCount64(admission.deferred) + " deferred; " +
        formatCount64(admission.published) + " published to the mid-quality output instead of counted\n";
    if (prefixMode_ || defersNonPassing()) {
        uint64_t recovered = 0;
        for (uint64_t count : recoveredCounts) { recovered += count; }
        report += string(prefixMode_ ? "Prefix matching: " : "Exact-key matching: ") + formatCount64(recovered) +
            " non-passing reads recovered into passing groups (deferred reads without a passing group are dropped)\n";
        if (admission.deferredDropped > 0) {
            report += string("Deferred table full (") + (retainsVariants() ?
                formatCount64(kDeferredEvidenceBytes >> 20) + " MiB of retained reads" : formatCount64(kDeferredKeyCap) + " keys") + "): " +
                formatCount64(admission.deferredDropped) + " non-passing reads with new keys not kept (counted as unmatched)\n";
        }
    }
    if (coarse_) {
        CoarseDereplication::Stats stats;
        size_t activePartitions = 0, searchGroups = 0;
        uint64_t largestPartition = 0;
        for (const auto& shard : shards_) {
            searchGroups += shard->coarse->groups().size();
            const auto& local = shard->coarse->stats();
            stats.queries += local.queries; stats.comparisons += local.comparisons;
            stats.fallbackQueries += local.fallbackQueries; stats.indexedSeeds += local.indexedSeeds;
            stats.exactHits += local.exactHits; stats.rejectedExactMisses += local.rejectedExactMisses;
            activePartitions += local.queries != 0;
            largestPartition = std::max(largestPartition, local.queries);
        }
        report += "Coarse search groups: " + formatCount64(searchGroups) + "; output " +
            (coarse_->storageOnly() ? "exact search-key dereplicates\n" : coarse_->outputsExact() ? "exact R1 dereplicates\n" : "coarse clusters\n");
        if (shards_.size() > 1) {
            report += "Coarse partitions: " + to_string(shards_.size()) + "; active " + to_string(activePartitions) +
                "; largest " + formatCount64(largestPartition) + " of " + formatCount64(stats.queries) + " queries\n";
        }
        report += "Coarse matching: " + formatCount64(stats.comparisons) + " bounded comparisons; " +
            formatCount64(stats.fallbackQueries) + " queries requiring seed fallback\n";
        report += "Coarse exact reuse: " + formatCount64(stats.exactHits) + " of " + formatCount64(stats.queries) +
            " queries; " + formatCount64(stats.queries - stats.exactHits - stats.rejectedExactMisses) + " candidate searches\n";
        report += "Coarse admission shortcut: " + formatCount64(stats.rejectedExactMisses) + " non-passing exact-key misses rejected without candidate search\n";
        if (coarseWorkPool_) {
            const auto work = coarseWorkPool_->stats();
            report += "Coarse work pool: " + formatCount64(work.batches) + " batches; " + formatCount64(work.tasks) +
                " partition tasks; " + formatCount64(work.helpedTasks) + " helped by another worker; peak active " +
                to_string(work.peakActive) + " partitions\n";
        }
        if (coarse_->reassigns()) {
            report += "Final reassignment: " + formatCount64(refinement.movedReads) + " reads in " +
                formatCount64(refinement.movedVariants) + " exact subclusters moved; " +
                formatCount64(refinement.searches.comparisons) + " representative comparisons\n";
        }
        if (coarse_->storesDiffs()) { report += "Nucleotide differences: " + diffPath + "\n"; }
        if (writeSubclusters_) { report += "Exact subcluster FASTA: " + subclusterPath + (pairedVariants ? ", " + subclusterMatePath : "") + "\n"; }
        if (coarse_->storesQualities()) {
            report += "Dereplicated HQ FASTQ (exact variants, selected qualities): " + outHQf + (pairedVariants ? ", " + outHQf_p2 : "") + "\n";
        }
    }
	report += "Search sequence length: avg " + avgSrchSeqLenStream.str() + " (" + formatCount64(srchSeqTotalLen) + " total over " + formatCount64(srchSeqCount) + " reads)\n";
	
		if (passed_hits > 0) {
	   report += N_notPassed + "/" + N_total + " not passing derep conditions (" + formatCount64(static_cast<uint64_t>(notPassedSize)) + " counts; "+ minCopiesStr;
		if (dereplicate_sample_specific) { report += " & sample specific restrictions"; }
	report += ")";
	}
		report += "\n"; 


		
		
	cerr  << report  << endl;


	if (this->b_derepPerSR) {
		return report;
	}

	int uneqCnts(0);
	//check counts_per_sample vector
	vector<string> SampleID = mf->SampleID;
	for (unsigned int i = 0; i<SampleID.size(); i++) {
		size_t j(0); bool detected = false;
		for (; j < barcode_number_to_sample_id_.size(); j++) {
			if (SampleID[i] == barcode_number_to_sample_id_[j]) { detected = true;  break; }
		}
        if (detected && mf->collectStatistics[0]->BarcodeDetected[i] != -1 &&
            static_cast<uint64_t>(mf->collectStatistics[0]->BarcodeDetected[i]) + recoveredCounts[j] != counts_per_sample[j]) {
			//cerr << "ERROR: Unequal counts for " << SampleID[i] << ": " << mf->collectStatistics.BarcodeDetected[i] << " vs " << counts_per_sample[j] << endl;
			uneqCnts++;
			
		} else if (!detected) {
            // j is past-the-end here; print the missing SampleID instead of indexing barcode_number_to_sample_id_
			cerr << "Could not detect sample_id_ " << SampleID[i] << " in ref set (check that mapping file is correctly formatted?).\n"; exit(87);
		}
	}
#ifdef DEBUG
	cerr << "Derep Fin" << endl;
#endif

	if (uneqCnts>0) {
		cerr << "Unequal counts in " << uneqCnts << " cases. \n";
		//exit(66);
	}
	return report;
}


bool Dereplicate::pass_deprep_conditions(shared_ptr<DNAunique> d) {
	if (minCopiesFastFailThreshold > 0 && d->totalSum() < static_cast<uint64_t>(minCopiesFastFailThreshold)) {
		return false;
	}
	size_t topSamples = 0;
	for (const DerepCopyRule& rule : minCopyRules) {
		topSamples = std::max(topSamples, static_cast<size_t>(rule.samples));
	}
	// Per-sample counts; the largest topSamples come first, in descending order.
	const vector<uint64_t> x = d->getDerepMapSort(topSamples);
	uint64_t sampleSum = 0;
	for (uint64_t count : x) { sampleSum += count; }
	if (x.empty()) { sampleSum = d->totalSum(); }
	for (const DerepCopyRule& rule : minCopyRules) {
		if (rule.copies <= 0) { return true; }
		if (rule.samples == 0) {// total across all samples
			if (sampleSum >= static_cast<uint64_t>(rule.copies)) { return true; }
		}
		// X:Y -> the Y-th most abundant sample still has at least X copies
		else if (x.size() >= static_cast<size_t>(rule.samples) &&
			x[rule.samples - 1] >= static_cast<uint64_t>(rule.copies)) {
			return true;
		}
	}
	return false;
}

void writeLog(string& logf) {
	ofstream of;
	of.open(logf.c_str());
	of.close();
}

string additionalFileName(const string& in){
	if (in.find(",") == std::string::npos){
		return additionalFileName2(in);
	}
	else {
		vector<string> vi = splitByCommas(in);
		string out = additionalFileName2(vi[0]);
		for (uint i = 1; i < vi.size(); i++){
			out += "," + additionalFileName2(vi[i]);
		}
		return out;
	}
}
string additionalFileName2(const string& in){
	size_t point = in.find_last_of('.');
	if (point == (size_t)-1){ return in + ".add"; }
	return in.substr(0, point) + ".add." + in.substr(point + 1);
}

string subfile(string x, string y) {
	if (x.empty() || y == "") { return x; }
	size_t pos = x.find_last_of(".");
	string s = "";
	if (pos != string::npos) {
		s= x.substr(0, pos) + "." + y + x.substr(pos);
	}
	else {
		s = x + "." + y;
	}
	return s;
}





////////////  UCLINKS  ///////////////
UClinks::~UClinks(){
	ucf.close();
	mapdere.close();
	if (merger != nullptr) {
		delete merger;
		merger = nullptr;
	}
}
UClinks::UClinks( OptContainer* cmdArgs):
	CurSetPair(-1),//maxOldDNAvec(20000),
	derepMapFile(""),
	bestDNA(0, NULL), oriKey(0), 
	mapLines(0),
	clusCnt(0), uclines(0),
	SEP(""), 
	UCread(false), pairsMerge(false), MAPread(false),
	b_derepAvailable(false),
	UPARSE8up(false), UPARSE9up(false), UPARSE11up(false),
	UpUcFnd(false), cdhit(false), vsearch(false),
	ucispaf(false),
	otuTerm("OTU"), otuOUTterm("OTU_"),
	RefDBmode(false), RefDBotuStart(-1),
	SeedsAreWritten(false),
	OTUmat(0), unregistered_samples(false),
	doChimeraCnt(false), OTUnumFixed(true),
	totalDerepCnt(0),
	qCovThr(0.8f), perIDmatch(97.f),
	b_merge_pairs_optiSeed_(false),merger(nullptr)

{
    const string subclusters = cmdArgs->get("-seedSubclusters", "0");
    if (subclusters != "0" && subclusters != "1") {
        throw invalid_argument("-seedSubclusters must be 0 or 1");
    }
    seedSubclusters = subclusters == "1";
    if (seedSubclusters && cmdArgs->get("-derep_map").empty()) {
        throw invalid_argument("-seedSubclusters requires the parent -derep_map");
    }
    if (seedSubclusters && !cmdArgs->get("-onlyPair").empty()) {
        throw invalid_argument("-seedSubclusters requires complete fragments; disable -onlyPair");
    }
	//read in UC file and assign clusters
	SEP = (*cmdArgs)["-sample_sep"];
	string str = (*cmdArgs)["-optimalRead2Cluster"];
	

	if (str.find(".paf") != string::npos) {
		ucispaf = true; UpUcFnd = true;
	}

	ucf.open(str.c_str(),ios::in);
	if (!ucf){
		UCread = true;
		cerr<<"Could not open uc file\n"<< str<<endl; exit(46);
	}
	if (cmdArgs->find("-derep_map") != cmdArgs->end()) {
		derepMapFile = (*cmdArgs)["-derep_map"];
	}
	if (cmdArgs->find("-minQueryCov") != cmdArgs->end()) {
		qCovThr = (float) atof((*cmdArgs)["-minQueryCov"].c_str());
	}
	if (cmdArgs->find("-id") != cmdArgs->end()) {
		perIDmatch =(float) atof((*cmdArgs)["-id"].c_str());
	}
	if (cmdArgs->find("-count_chimeras") != cmdArgs->end() &&
		(*cmdArgs)["-count_chimeras"] == "T") {
		doChimeraCnt = true;
	}
	if (cmdArgs->find("-merge_pairs_seed") != cmdArgs->end() &&
		(*cmdArgs)["-merge_pairs_seed"] == "1") {
		b_merge_pairs_optiSeed_ = true;
	}

	//set version of mapping/clustering
	if ((*cmdArgs)["-uparseVer"] != "") {
		if ((*cmdArgs)["-uparseVer"] == "N11") {//UNOISE v11
			UPARSE8up = true;
			otuTerm = "Zot";
			otuOUTterm = "Zotu";
		}
		else if ((*cmdArgs)["-uparseVer"] == "cdhit") {
			cdhit = true; UpUcFnd = true; UPARSE8up = true;
		} else if ((*cmdArgs)["-uparseVer"] == "vsearch") {
			UpUcFnd = true; vsearch = true;
		} else {
			int upVer = atoi((*cmdArgs)["-uparseVer"].c_str());
			if (upVer >= 8 && upVer < 9) {
				UPARSE8up = true; UpUcFnd = true;
			}
			else if (upVer >= 9 && upVer < 11) {
				UPARSE8up = true; UPARSE9up = true; UpUcFnd = true;
			}
			else if (upVer == 9966) {//cdhit code
				cdhit = true; UpUcFnd = true;
			}
			else if (upVer >= 11) {
				UPARSE8up = true; UPARSE9up = true; UPARSE11up = true; UpUcFnd = true;
				otuTerm = "otu";
			}
			
		}
	}
}


void UClinks::readDerepInfo(const string dereM) {
	mapdere.open(dereM.c_str(), ios::in);
	if (!mapdere) {
		cerr << "Can't open " << dereM << ". \nAborting\n"; exit(54);
	}
	int cnt(0);
	totalDerepCnt = 0;
	b_derepAvailable = true;
	string line("");
	//only read header
	while (!safeGetline(mapdere, line).eof()) {
		//	while(getline(in,line,'\n')) {
		if (line.length() < 3) { continue; }
		string segments;
		stringstream ss;
		ss << line;
		int tbcnt = -1;

		//(*cmdArgs)["-i_MID_fastq"]
		while (getline(ss, segments, '\t')) {
			tbcnt++;
			if (cnt == 0) { //search for header
				if (tbcnt == 0) {
					if (segments != "#SMPLS") { cerr << "First line in dereplicate map has to start with #SMPLS, corrupted file.\n" << dereM << endl; exit(98); }
					continue;
				}
				vector<string> spl = header_string_split(segments, ":");
				int ccnt = stoi(spl[0]);
				SmplIDs[spl[1]] = ccnt;				
			} else {//actual derep info per sorted line
				cerr << "wrong mapping reading\n"; exit(64);
			}
		}
		if (cnt == 0) {
			OTUmat.resize(SmplIDs.size(), vector<matrixUnit>(clusCnt + 1, 0));
		}
		break;
	}
	cout << "Found " << SmplIDs.size() << " samples in derep.map\n";
}

// Shared map reader: HQ input consumes one row at a time; reconstructed input
// loads the parent counts once and joins them to mappings independently of order.
bool UClinks::nextDerepLine(string& line) {
    if (MAPread) { return false; }
    while (safeGetline(mapdere, line)) {
        if (!line.empty() && line[0] != '#') { return true; }
    }
    mapdere.close();
    MAPread = true;
    return false;
}

// Abundance annotation of a dereplicate header: ";size=N;" (default) or a
// trailing "_N" (-dere_size_fmt 1, which seed extension used to reject). 0 if none.
static uint64_t dereplicateHeaderAbundance(const string& header) {
    const size_t size = header.find(";size=");
    if (size != string::npos) { return stoull(header.substr(size + 6)); }
    const size_t underscore = header.find_last_of('_');
    if (underscore != string::npos && underscore + 1 < header.size()) {
        const string digits = header.substr(underscore + 1);
        if (is_digits(digits)) { return stoull(digits); }
    }
    return 0;
}

uint64_t UClinks::parseDerepLine(const string& line, shared_ptr<DNAunique> d, bool strictCounts) {
    const size_t tab = line.find('\t');
    const string header = line.substr(0, tab);
    if (tab == string::npos) { throw runtime_error("Missing sample counts in derep.map: " + line); }
    const uint64_t declared = dereplicateHeaderAbundance(header);
    uint64_t count = 0;
    string entry;
    stringstream samples(line.substr(tab + 1));
    while (getline(samples, entry, '\t')) {
        if (entry.empty()) { continue; }
        const size_t colon = entry.find(':');
        if (colon == string::npos) { throw runtime_error("Invalid sample count in derep.map: " + line); }
        const int sample = stoi(entry.substr(0, colon));
        const uint64_t copies = stoull(entry.substr(colon + 1));
        if (copies > numeric_limits<uint64_t>::max() - count) {
            throw overflow_error("derep.map abundance overflow: " + header);
        }
        d->incrementSampleCounterBy(sample, copies);
        count += copies;
    }
    if (count != declared) {
        const string error = "Mapping file unique abundance reconstruction: failed to find correct count (" +
            to_string(count) + " vs " + to_string(declared) + "):\n" + line;
        if (strictCounts) { throw runtime_error(error); }
        cerr << "ERROR: " << error << endl;
    }
    totalDerepCnt += count;
    return count;
}

uint64_t UClinks::oneDerepLine(shared_ptr<DNAunique> d) {
    string line;
    if (!nextDerepLine(line)) { return 0; }
    const string header = line.substr(0, line.find('\t'));
    if (!d->sameHead(header)) {
        cerr << header << " is not " << d->getId() << endl;
        exit(738);
    }
    return parseDerepLine(line, d, false);
}

string UClinks::seedParentKey(string id) {
    id = id.substr(0, id.find_first_of(" \t\r\n"));
    removeSizeStr(id);
    if (!SEP.empty()) { removeSampleID(id, SEP); }
    return id;
}

// Count-only parents remain in unusedDNA until a primary or additional mapping
// consumes them, retaining the existing .rest/reference count-recovery behavior.
void UClinks::finishMAPfile() {
    if (MAPread) {
        if (seedSubclusters && subclusterMappingsLoaded) {
            for (const auto& [id, parent] : subclusterParents) {
                if (parent.lastVariant && parent.observed != parent.expected) {
                    throw runtime_error("Incomplete subclusters for parent " + id + ": observed " +
                        to_string(parent.observed) + " fragments, expected " + to_string(parent.expected));
                }
            }
        }
        return;
    }
    string line;
    uint64_t finishMapCnt = 0;
    while (nextDerepLine(line)) {
        const string header = line.substr(0, line.find('\t'));
        auto d = make_shared<DNAunique>("", header);
        d->seal();
        const uint64_t count = parseDerepLine(line, d, true);
        finishMapCnt += count;
        const string key = seedSubclusters ? seedParentKey(header) : header.substr(0, header.find(' '));
        if (seedSubclusters) {
            if (!subclusterParents.emplace(key, SubclusterParent{count}).second) {
                throw runtime_error("Duplicate parent in derep.map: " + key);
            }
        }
        unusedDNA[key] = d;
    }
    cerr << "Found " << totalDerepCnt << " counts in derep.map, " << finishMapCnt << " counts secondary\n";
}

shared_ptr<DNAunique> UClinks::prepareSeedCandidate(shared_ptr<DNA> first,
        shared_ptr<DNA> second, Filters* fil) {
    // HQ and reconstructed reads already have their physical preprocessing cuts.
    auto candidate = make_shared<DNAunique>(first, -1);
    if (second) {
        candidate->attachPair(make_shared<DerepRead>(second));
        if (merger) {
            ++seedMergeAttempts_;
            if (merger->findSeedForMerge(candidate, second)) { ++seedMergeSuccesses_; }
        }
        if (merger && candidate->merge_seed_pos_ >= 0) {
            auto merged = merger->merge(candidate, second);
            if (merged && fil->checkMergedDNA(merged, true)) {
                candidate->attachMerge(make_shared<DerepRead>(merged));
            }
        }
    }
    return candidate;
}

void UClinks::findSubclusterSeeds(shared_ptr<InputStreamer> input, bool readFQ, Filters* fil) {
    if (!readFQ || !input->qualityPresent() || input->hasMIDseqs()) {
        throw invalid_argument("-seedSubclusters requires reconstructed single or paired FASTQ with qualities");
    }
    if (!b_derepAvailable) { throw runtime_error("-seedSubclusters requires initialized parent counts and fallback seeds"); }
    if (!subclusterMappingsLoaded) {
        finishMAPfile();
        string query, target;
        float identity = 0.f;
        vector<int> clusters;
        while (getMAPPERline(query, target, identity, clusters, false)) {
            if (clusters.empty()) { continue; }
            auto parent = subclusterParents.find(seedParentKey(query));
            int counted = 0;
            if (parent != subclusterParents.end() && uclInOldDNA_simple(query, clusters, counted)) {
                parent->second.clusters = clusters;
                parent->second.identity = identity;
            }
        }
        subclusterMappingsLoaded = true;
    }
    // Only assignments/count metadata are retained per parent. Each full exact
    // variant competes in export order, with its linked mate and quality vectors.
    while (true) {
        const auto reads = input->getDNApairs();
        if (!reads[0]) { break; }
        const string id = reads[0]->getShortId();
        if (reads[1] && id != reads[1]->getShortId()) {
            throw runtime_error("Subcluster mate IDs and abundances must match: " + id + " versus " + reads[1]->getShortId());
        }
        const size_t size = id.find(";size=");
        const string label = id.substr(0, size);
        const size_t suffix = label.rfind(".sub");
        const string number = suffix == string::npos ? "" : label.substr(suffix + 4);
        const string copies = size == string::npos ? "" : id.substr(size + 6, id.find(';', size + 6) - (size + 6));
        if (suffix == string::npos || suffix == 0 || number.empty() || copies.empty() ||
                !is_digits(number) || !is_digits(copies)) {
            throw runtime_error("Invalid reconstructed subcluster ID (expected parent.subN;size=N;): " + id);
        }
        const uint64_t variant = stoull(number), count = stoull(copies);
        auto entry = subclusterParents.find(seedParentKey(label.substr(0, suffix)));
        if (entry == subclusterParents.end()) { throw runtime_error("No derep.map parent for subcluster " + id); }
        auto& parent = entry->second;
        if (!count || variant != parent.lastVariant + 1 || count > parent.expected - parent.observed) {
            throw runtime_error("Invalid subcluster order or abundance for " + id);
        }
        parent.lastVariant = variant;
        parent.observed += count;
        if (parent.clusters.empty()) { continue; }
        auto candidate = prepareSeedCandidate(reads[0], reads[1], fil);
        candidate->setTempFloat(parent.identity);
        besterDNA(parent.clusters, candidate, fil, false, count);
    }
    cerr << "Associated " << totalDerepCnt << ", " << OTUmatSum() << " counts to primary clusters\n";
}

//functions selects optimal read to represent OTU
void UClinks::findSeq2UCinstruction(shared_ptr<InputStreamer> IS, bool readFQ, 
	Filters* fil){
    if (seedSubclusters) {
        findSubclusterSeeds(IS, readFQ, fil);
        return;
    }
	if (UCread){return;}

	if (IS->hasMIDseqs()){
		CurSetPair = 0;
	} else {
		CurSetPair = -1;
	}
#ifdef DEBUG
	cerr << "UC seed extension" << endl;
#endif
	uint64_t UCcnt = 0;
	std::list<string>::iterator listIT;
	shared_ptr<DNAunique> match(nullptr); shared_ptr<DNA>  match2(nullptr);
	string segs("");
	string segs2("");
	float perID(0.f);
	vector<int> curCLID(0,0);
	//bool sync(false); // syncing of 2 read pairs; not implemented for this function yet


	//goes through UC file
	while ( getMAPPERline(segs, segs2, perID, curCLID, !b_derepAvailable) ) {
		//int subcnt = 0;
		if ( curCLID.size() == 0 ) { continue; }
		//first check if we already have this read in oldDNA (e.g. from previous lines), if so, 1) integrate 2) skip to next line
		if ( uclInOldDNA(segs, curCLID, perID, fil) ) {
			curCLID.resize(0);
			continue;
		}
        // for each mapping entry, allow reading DNA records until a match is found or EOF
		bool cont_local = true;
		while(cont_local){//match to DNA
			shared_ptr<DNA> tmpDNA = IS->getDNA(0);
			if (tmpDNA == NULL) { cont_local = false; break; }//signal that at end of file
			match2 = IS->getDNA(1);

			match = prepareSeedCandidate(tmpDNA, match2, fil);

			//assummes in original implementation, that we can get derep.map lines with the same 
			//ordering as fq derep (which works normally, just not for dada2 mode)
			UCcnt+= oneDerepLine(match);
			string curID("");
            if (match != nullptr) {
                curID = match->getId(); 
                curID = curID.substr(0, curID.find_first_of(' ')); 
                match->setTempFloat(perID); 
            } 
            //check if dnaTemp1 a) matches id_ b) is better 
            if (curID != segs) { 
                if (match != nullptr) {
					//block to store unused DNA & find this id_ in this block
                    unusedDNA[curID] = match;
					//oldDNA[DNAunusedPos] = match;	//oldDNA2[DNAunusedPos] = match2;
                    //DNAunusedPos++;
				}
            } else {//find if current DNA rep needs to be replaced with a better read (pair) 
				//assign % identity score to DNA object
				cdbg( "UC Hit");				
				besterDNA(curCLID, match, fil);
				curCLID.resize(0);
				break;
			}
			
			//if (tdn!=NULL && ch1 != tdn->isGreenQual()){cerr<<"isGreenQual is != ch1! Aborting..\n";exit(12);}
		}
	}
#ifdef DEBUG
	cerr << "UC seed initialized" << endl;
#endif
	
	
	//cerr << "\nRead " << totalDerepCnt << " read counts from derep.map\n";
	cerr << "Associated " << totalDerepCnt << ", "<< OTUmatSum() <<" counts to primary clusters\n";

}

void UClinks::finishUCfile(Filters* fil, string addUC, bool bSmplHd){
	if (UCread){
		addUCdo(addUC,bSmplHd);
		UCread = true;
		return;
	}
	string segs;
	string segs2;
	float perID(0.f);
	vector<int> curCLID(0);

	while ( getMAPPERline(segs, segs2, perID, curCLID, !b_derepAvailable) ) {
		if ( uclInOldDNA(segs, curCLID, perID, fil) ) {
			curCLID.resize(0);
			continue;
		}
		cerr << segs << " ";
	}
	
	addUCdo(addUC, bSmplHd);
	UCread = true;
}
void UClinks::addUCdo(string addUC,bool SmplHd) {
	if (addUC == "") {
		return;
	}
	string segs;
	string segs2;
	float perID;
	int cntsAddUC=0;
	vector<int> curCLID(0);
	ucispaf = false;
	if (addUC.find(".paf") != string::npos) {
		ucispaf = true; UpUcFnd = true;
	}
	ucf.clear();
	ucf.open(addUC.c_str(), ios::in);
	if (!ucf) {
		// Nothing to add; parsing a failed stream could not end.
		cerr << "Could not find additional uc file: " << addUC << endl;
		ucf.clear();
		UCread = true;
		return;
	}
	std::cerr << "Reading " << addUC << endl;
	UCread = false;
	if (SmplHd){
		while (getMAPPERline(segs, segs2, perID, curCLID, SmplHd)) { 
			if (curCLID.size() > 0) { cntsAddUC++; }
			curCLID.resize(0); 
		}
		
	}else {//complicated..
		int cnt(0);
		while (getMAPPERline(segs, segs2, perID, curCLID, SmplHd)) { 
			//find segs in remaining dereps
			cnt++; //if (cnt < 61403) { continue; }
			if (curCLID.size() == 0) { continue; }
			if (uclInOldDNA_simple(segs, curCLID, cntsAddUC)) {
				curCLID.resize(0);
				continue;
			}

			curCLID.resize(0); 
		}
	}

	cout << "Added " << cntsAddUC << " counts from uc add file\n";

}

//used for specific situation where only empty DNA string was read with derep info attached
bool UClinks::uclInOldDNA_simple(const string& segs,const vector<int>& curCLID, int& countsAdd) {
	//check if DNA is in group of unmatched DNA's ?
	if (segs == ""){ return false; }//empty ucl line
	auto unusedIT = unusedDNA.find(seedSubclusters ? seedParentKey(segs) : segs);
	if (unusedIT != unusedDNA.end()){// found something
		//int mID = (*unusedIT).second;
		//if (mID >= (int)oldDNA.size()) { cerr << "SEC MID too high\n"; exit(55); }
		if (unusedIT->second == nullptr){ return false; }
		//give sequence a chance to be selected
        addSeedCounts(curCLID, unusedIT->second);
        countsAdd += (int)curCLID.size();
		//was matched once to an OTU seed. Even if several later matches, doesn't mater - delete

		unusedDNA.erase(unusedIT);
//		delete oldDNA[mID]; if (oldDNA2[mID] != NULL){ delete oldDNA2[mID]; }

		//oldDNA.erase(mID); //oldDNA2.erase(mID);// = NULL;oldDNA2[mID] = NULL;
		return true;
	}
	return false;
}
bool UClinks::uclInOldDNA(const string& segs,const vector<int>& curCLID, float perID,
	Filters* fil) {
	//check if DNA is in group of unmatched DNA's ?
	if (segs == ""){ return false; }//empty ucl line
	auto unusedIT = unusedDNA.find(seedSubclusters ? seedParentKey(segs) : segs);
	if (unusedIT != unusedDNA.end()){// found something
		//int mID = (*unusedIT).second;
		if (unusedIT->second == nullptr) { return false; }
		//give sequence a chance to be selected
		unusedIT->second->setTempFloat(perID);
		besterDNA(curCLID, unusedIT->second, fil);
		//remove all trace
		unusedDNA.erase(unusedIT);
		return true;
	}
	return false;
}

void UClinks::resetMarks() {
	if (cdhit || vsearch) {
		UpUcFnd = false;
		cdhit = false; vsearch = false;
	}
}

bool UClinks::getTMPmapperLine(string& line) {
	if (mapLines.size() > 0) {
		line = mapLines.front();
		mapLines.pop_front();
	}
	else {
		if (!readMapperFileLine(line)) {
			return false;
		}
	}
	return true;
}

// One line of the mapper file without a trailing '\r' (CRLF files made the
// cd-hit representative test `line.back() == '*'` fail forever).
bool UClinks::readMapperFileLine(string& line) {
	if (!getline(ucf, line, '\n')) { return false; }
	if (!line.empty() && line.back() == '\r') { line.pop_back(); }
	return true;
}

// Normalizes an OTU/target name the same way for lookups and for the
// perfect-chimera set: no description, size, sequence or centroid tags, no
// sample prefix.
string UClinks::clusterKey(string name) {
	name = name.substr(0, name.find_first_of(' '));
	removeSizeStr(name);
	removeSeqStr(name);
	removeCentrStr(name);
	removeSampleID(name, SEP);
	return name;
}

bool UClinks::getMAPPERline(string& segs, string& segs2,float& perID, 
	vector<int>& curCLID,  bool addFromHDstring) {
	curCLID.clear();
	//converts derep reads mapping to OTUs to counts
	//seg: read ID, segs2: OTU repr read, perID: match to repr read


	//close all file streams
	if (UCread) { resetMarks(); UpUcFnd = false; return false; }
	// Buffered cd-hit members are still pending when the file hit EOF on the
	// representative line; finish them before closing.
	if (ucf.eof() && mapLines.empty()){
		resetMarks();
		UCread=true;		ucf.close();		return false;
	}
	float qCov = 1.f;
	string line; 
	std::unordered_map<string, int>::iterator itCL;

	//check if any lines need to be backworked
	while (this->getTMPmapperLine(line)) {
	//	cerr<<line<<endl;
		uclines++;
		if (line.length() <= 1){ continue; }	
      if (!UpUcFnd) {
			// reset
			vsearch = false; cdhit = false; UPARSE8up = false; UPARSE9up = false; UPARSE11up = false;
			if (line[0] == '>') {
				cdhit = true;
			} else if (line.substr(1, 1) == "\t"){ UPARSE8up = false; 
			}else {	UPARSE8up = true;	}

			if (UPARSE8up){//could also be uparse9

				stringstream ss2;
				ss2 << line;   int tabCnt = 1;
				while (getline(ss2, segs, '\t')) {
					tabCnt++;
				}
				if (tabCnt >= 3) {
					UPARSE9up = true;
					cerr << "Switching to Uparse 10+ style map file.\n";
				}
			}
			UpUcFnd = true;
		}
		bool chimera = false;
		vector<string>tarsV; //saves hits to OTUs

		if (cdhit) {
			// Member name between "nt, >" and "..."; a malformed line used to
			// yield garbage names (or out_of_range) from unchecked positions.
			auto memberName = [](const string& clstr, size_t& end) {
				const size_t pos = clstr.find("nt, >");
				end = pos == string::npos ? string::npos : clstr.find("...", pos + 5);
				if (end == string::npos) { throw runtime_error("Malformed CD-HIT .clstr line: " + clstr); }
				return clstr.substr(pos + 5, end - pos - 5);
			};
			if (line[0] == '>') {
				//get next line with actual mapping/clustering
				if (!readMapperFileLine(line) || line.empty()) {
					throw runtime_error("CD-HIT .clstr cluster header without members at end of file");
				}
				repFound = false;
				segs2 = "";
			}

			//meachnism has to be slightly different for cdhit, have to save segs2 between lines
			if (!repFound  ){//&& line.back() == '*') {//report representative gene
				bool isRep = !line.empty() && line.back() == '*';
				while (!isRep) {
					mapLines.push_back(line);
					// getline at EOF left `line` unchanged: the loop never ended.
					if (!readMapperFileLine(line)) {
						throw runtime_error("CD-HIT .clstr cluster without representative (*) at end of file");
					}
					isRep = !line.empty() && line.back() == '*';
					if (!line.empty() && line[0] == '>') {//should never get here
						cerr << "CD-HIT .clstr contained seq cluster with no rep (prev lines):\n" << line << endl;
						exit(163);
					}
				}
				size_t end = 0;
				segs = memberName(line, end);
				repFound = true;
				segs2 = segs;
				perID = 100.f;
				
			} else {
				size_t end = 0;
				segs = memberName(line, end);
				//find ID
				const size_t pos = line.find("at +/", end + 3);
				if (pos == string::npos || line.size() < pos + 6) { throw runtime_error("Malformed CD-HIT .clstr identity: " + line); }
				string ID = line.substr(pos + 5, line.length() - pos - 6);
				perID = stof(ID);
			}
		} else {

			stringstream ss;
			ss << line;
			//2 ways to get to a) hit info b) query & otu
			if (ucispaf) {//minimap2 .paf file..
				getline(ss, segs, '\t');
				string query = segs;
				getline(ss, segs, '\t');
				float qLen = (float)atof(segs.c_str());
				getline(ss, segs, '\t');
				const float qStart = (float)atof(segs.c_str());
				getline(ss, segs, '\t');
				const float qEnd = (float)atof(segs.c_str());
				getline(ss, segs, '\t');//strand
				getline(ss, segs2, '\t');//tarname
				getline(ss, segs, '\t');
				float tLen = (float) atof(segs.c_str());
				for (uint i = 0; i < 2; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				getline(ss, segs, '\t');
				float matches = (float) atof(segs.c_str());
				getline(ss, segs, '\t');
				float mapLen = (float) atof(segs.c_str());
				perID = matches / mapLen *100.f;

				// Query coverage is the aligned query span, not the alignment block
				// length (which includes target-relative gaps and can exceed the query).
				qCov = qLen > 0.f ? (qEnd - qStart) / qLen : 0.f;
				// Secondary (tp:A:S) and inversion (tp:A:I) records repeat a query;
				// counting them would count one read several times.
				const bool secondary = line.find("\ttp:A:S") != string::npos || line.find("\ttp:A:I") != string::npos;
				if (secondary || perID < perIDmatch || qCov < qCovThr) {
					//reject this hit
					return true;
				}
				segs = query;
			} else if (vsearch) {
				bool isSeed = line.substr(0, 1) == "S";
				if (line.substr(0, 1) != "H" && !isSeed) {
					continue;
				}
				for (uint i = 0; i < 4; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				perID = (float)atof(segs.c_str());
				for (uint i = 0; i < 5; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				getline(ss, segs2, '\t');
				if (isSeed) { segs2 = segs; perID = 100.f; }

			} else if (!UPARSE8up) { //uparse 7
				if (line.substr(0, 1) != "H" ) {
					continue;
				}
				for (uint i = 0; i < 4; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				perID = (float)atof(segs.c_str());
				for (uint i = 0; i < 5; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				getline(ss, segs2, '\t');
			}
			else if (!UPARSE9up) { // uparse 8
			 //query first entry
				string tmp;
				getline(ss, segs, '\t');//0
				getline(ss, tmp, '\t');//1
				//should be "match"
				if (tmp == "chimera") {
					if (!doChimeraCnt) { continue; }
					chimera = true;
				}// segs = ""; continue;}
				if (tmp == "otu" || tmp == "OTU") {
					segs2 = segs;
					perID = 100.f;
				}
				else {
					getline(ss, tmp, '\t');//2
					perID = (float)atof(tmp.c_str());
					//indicator if hit
					//OTU last entry
					getline(ss, segs2, '\t');//3
					getline(ss, segs2, '\t');//4
				}

			}
			else { //UP9, uparse 10, uparse 11, dada2 fake .uc
				  //query first entry
				string tmp;
				getline(ss, segs, '\t');//0
				getline(ss, tmp, '\t');//1
									   //should be "match"
				if (tmp == "chimera") {
					continue;
				}
				else if (tmp == "noisy_chimera" || tmp == "good_chimera") { //tmp == "perfect_chimera" ||
					if (!doChimeraCnt) { continue; }
					chimera = true;
//M04428:252:000000000-BM3M5:1:1109:4503:15901;size=52;   noisy_chimera   
// dqt=20;dqm=2;div=9.0;segs=2;parents=M04428:252:000000000-BM3M5:1:2110:24763:14436(1-78/0)+M04428:252:000000000-BM3M5:1:1116:5359:16203(78-186/0);

					getline(ss, tmp, '\t');
					// "parents=" is 8 characters; +4 pointed into the keyword and never matched an OTU.
					const size_t parents = tmp.find("parents=");
					if (parents == string::npos) { continue; }
					size_t p1(parents + 8);
					size_t p2(tmp.find(");", p1) + 2);
					segs2 = tmp.substr(p1, p2 - p1 - 1);
					perID = 100.f;
				}
				else if (tmp == "perfect_chimera") {
					// Same normalization as the target lookup below, so the
					// "unknown OTU" suppression actually matches.
					perfectChims.insert(clusterKey(segs));
					continue;
				}// segs = ""; continue;}
				if (tmp.substr(0, 3) == otuTerm) {
					segs2 = segs;
					perID = 100.f;
				}
				else if (!chimera) {//match case (a chimera line has no top= column; parsing it overwrote the parents)
					getline(ss, tmp, '\t');//2
					//dqt=1;top=GZV0ATA01ANJXZ;size=14;(99.6%);
					size_t p1(tmp.find("top=") + 4);//4
					size_t p2(tmp.find(";(", p1) + 2);
					size_t p3(tmp.find("%);", p2));
					segs2 = tmp.substr(p1, p2 - p1 - 1);
					//string xx = tmp.substr(p2, p3 - p2);
					perID = (float)atof(tmp.substr(p2, p3 - p2).c_str());
					if (false && chimera) {//just use up9 top hit
	//					p1(tmp.find(";top=") + 5,p2);
	//					p2(tmp.find(";(", p1) + 2);

					}
				}
			}
		}



		//convert mapping to otu information
		//remove spaces
		segs = segs.substr(0,segs.find_first_of(' '));
		//also remove sample identifier in string
		string smplID = "";
		removeSampleID(segs, SEP, smplID);

		//finished parsing, now reformat to get matching sample

		if (chimera){
			if (UPARSE9up) {
				tarsV = splitByComma(segs2, false, '+');
				for (uint kk = 0; kk < tarsV.size(); kk++) {
					tarsV[kk] = tarsV[kk].substr(0, tarsV[kk].find_last_of("("));
				}
			}
			else if (UPARSE8up) {
				tarsV = splitByComma(segs2, false, '+');
				for (uint kk = 0; kk < tarsV.size(); kk++) {
					tarsV[kk] = tarsV[kk].substr(0, tarsV[kk].find_last_of("("));
				}
			}
		} else {
			tarsV.push_back(segs2);
		}
		// -1 marks unknown targets; zero-filling silently assigned them to OTU 0.
		curCLID.assign(tarsV.size(), -1);

		for ( uint kk = 0; kk < tarsV.size(); kk++ ) {
			//goes over all chimeric hits
			string oriClKey (tarsV[kk]);
			segs2 = clusterKey(oriClKey);
			itCL = seq2CI.find(segs2);

			if ( itCL == seq2CI.end() ) {//new OTU?
				if (OTUnumFixed) {
					//cerr << "XX\n";
					if (perfectChims.find(segs2) == perfectChims.end()) {
						cerr << "Unkown OTU entry found:" << segs2 << endl;
					}
				}	else {
					//not found in known clusters.. create entry
					bestDNA.push_back(make_shared<betterSeedStruct>());
					//bestDNA2.push_back(NULL);
					oriKey.push_back(oriClKey);
					clusCnt = (int)bestDNA.size() - 1;
					curCLID[kk] = clusCnt;
					seq2CI[segs2] = clusCnt;
				}
			} else {//cluster exists
				curCLID[kk] = ((*itCL).second);
			}
		}
		curCLID.erase(std::remove(curCLID.begin(), curCLID.end(), -1), curCLID.end());
#ifdef matrix_sum
		// A chimera's count is split over its known targets, as addSeedCounts
		// does for dereplicated reads (it was divided by all targets here, so an
		// unknown parent silently lost that share).
		if (addFromHDstring && smplID != "") {
			const matrixUnit splChim = (matrixUnit)curCLID.size();
			for (int cluster : curCLID) { add2OTUmat(smplID, cluster, splChim); }
		}
#endif
		return true;
	}
	// The file ended or failed: finish like EOF. Returning true here let a
	// stream that failed without EOF spin the callers' while loops forever.
	resetMarks();
	UCread = true;		ucf.close();
	return false;
}
void UClinks::writeOTUmatrix(string outf) {
	cerr << "Writing OTU matrix with "<< SmplIDs.size() << " samples and " << seq2CI.size() << " OTUs to " << outf << endl;
	this->setOTUnms();
	ofstream MA;
	MA.open(outf);
	//first write all sample IDs
	std::unordered_map<string, int>::iterator OTUid;
	MA << "OTU";
	matrixUnit totalCnts = 0;
	
	vector<string> smpls(SmplIDs.size(),"");
	unordered_map<string, int>::iterator smplNit;
	for (smplNit = SmplIDs.begin(); smplNit != SmplIDs.end(); smplNit++) {
		smpls[smplNit->second] = smplNit->first;
	}
	
	for (size_t i = 0; i < smpls.size(); i++) {
		MA<<"\t"<<smpls[i];
	}
	//write OTU values
	for (OTUid = seq2CI.begin(); OTUid != seq2CI.end(); OTUid++) {
		MA << endl << oriKey[OTUid->second];
		int rowI = OTUid->second;
		for (size_t i = 0; i < smpls.size(); i++) {
			matrixUnit val = OTUmat[ SmplIDs[smpls[i]] ][rowI];
			totalCnts += val;
			if ( isnan(val) ) {
				MA << "\t0";
			} else {
				// Print as an integer: a double streams with 6 significant digits
				// (1234567 became 1.23457e+06 above one million reads).
				MA << "\t" << static_cast<unsigned long long>(std::llround(std::max(0.0, (double)val)));
			}
		}
	}
	MA << endl;
	MA.close();
	cerr << "Recruited " << static_cast<unsigned long long>(std::llround(std::max(0.0, (double)totalCnts))) << " reads in OTU matrix\n";
}
void UClinks::setOTUnms(){
	int cnt = 0;
	std::unordered_map<string, int>::iterator OTUid;
    for (OTUid = seq2CI.begin(); OTUid != seq2CI.end(); OTUid++) {
		string newlySetID = otuOUTterm + std::to_string(cnt);
		oriKey[OTUid->second] = newlySetID;
		cnt++;
	}

}
void UClinks::add2OTUmat(const string& smplID, int curCLID, matrixUnit spl) {

	std::unordered_map<string, int>::iterator smplNit = SmplIDs.find(smplID);
	if (smplNit == SmplIDs.end()) {//create entry & expand matrix
		SmplIDs[smplID] = (int) OTUmat.size();
		OTUmat.push_back(vector<matrixUnit>(clusCnt + 1, (matrixUnit)0));
		smplNit = SmplIDs.find(smplID);
		cerr << "New sample_id_ id_ in uc file detected, that is not present in map: \"" << smplID << "\"\n";// endl;
		unregistered_samples = true;
	}
	//easy, now add in
	if ( spl <1 ) { spl = 1; }
	OTUmat[(*smplNit).second][curCLID] += (matrixUnit)1 / spl;

}
void UClinks::add2OTUmat(shared_ptr<DNAunique> d, int curCLID, matrixUnit rep) {
	if (d == NULL) { cerr << " add2OTUmat::d is NULL\n"; exit(85); }
	const read_occ& map = d->getDerepMap();
	if ( rep <1 ) { rep = 1; }
	for (const auto& entry : map) {
		if (entry.second == 0 || entry.first < 0 || static_cast<size_t>(entry.first) >= OTUmat.size()) { continue; }
		OTUmat[entry.first][curCLID] += static_cast<matrixUnit>(entry.second) / rep;
	}
}

void UClinks::setupDefSeeds(shared_ptr<InputStreamer> FA, const vector<string>& smpls) {
	bool contRead = true; //bool sync(false);
	while (contRead) {
		shared_ptr<DNA> tmpDNA = FA->getDNA( 0);
		if (tmpDNA == NULL) { contRead = false; break; }
		tmpDNA->setAllQual(30);
		shared_ptr<DNAunique> tmp = make_shared<DNAunique>(tmpDNA, -1);
//		delete tmpDNA;
		//second pair_
		//shared_ptr<DNA> tmp2 = FA->getDNA(contRead, 1);
		string segs2 = tmp->getId();
		string oriClKey = segs2;
		//remove ;size= argument
		removeSizeStr(segs2);
		//also remove sample identifier in string
		removeSampleID(segs2, SEP);
		if (vsearch) {
			removeSeqStr(segs2);
			removeCentrStr(segs2);
		}

		//cluster should not exist, test
		std::unordered_map<string, int>::iterator itCL;
		itCL = seq2CI.find(segs2);

		if (itCL == seq2CI.end()) {
			//not found in known clusters.. create entry
			tmp->setTempFloat(100.f);
			shared_ptr< betterSeedStruct> tmpDNAst = make_shared<betterSeedStruct>(tmp);
			bestDNA.push_back(tmpDNAst);

			//bestDNA2.push_back(NULL);
			oriKey.push_back(oriClKey);
			clusCnt = (int)bestDNA.size() - 1;
			seq2CI[segs2] = clusCnt;
		} else {//cluster exists
			cerr << "Found double id_ " << segs2 << endl;//
			if (false) {
				cerr << "Aborting.." << endl;
				exit(74);
			}
		}
	}
	//std::map<string, int>::iterator smplNit;
	if (derepMapFile != "") {
		readDerepInfo(derepMapFile);
	} else {
//		const vector<string>& smpls = fil->SampleID;
		for (uint i = 0; i < smpls.size(); i++) {
			SmplIDs[smpls[i]] = (int)OTUmat.size();
			if (i != OTUmat.size()) { cerr << "Err in setupDefSeeds\n"; exit(453); }
			OTUmat.push_back(vector<matrixUnit>(clusCnt + 1, (matrixUnit)0));
		}
	}
}

void UClinks::addDefSeeds(shared_ptr<InputStreamer> FA, Filters* fil) {
	bool contRead = true;
	int addCnt = 0; //bool sync(false);
	while (contRead) {
		shared_ptr<DNA> tmpDNA = FA->getDNA( 0);
		if (tmpDNA == NULL) { contRead = false; break; }
		shared_ptr<DNAunique> tmp = make_shared<DNAunique>(tmpDNA, -1);
		//delete tmpDNA;
		//second pair_
		//shared_ptr<DNA> tmp2 = FA->getDNA(contRead, 1);
		string segs2 = tmp->getId();
		string oriClKey = segs2;
		//remove ;size= argument
		removeSizeStr(segs2);
		removeSampleID(segs2, SEP);
		if (vsearch) {
			removeCentrStr(segs2);
			removeSeqStr(segs2);
		}

		//cluster should not exist, test
		std::unordered_map<string, int>::iterator itCL;
		itCL = seq2CI.find(segs2);

		if (itCL == seq2CI.end()) {
			//not found in known clusters.. create entry
			tmp->setTempFloat(100.f);
			shared_ptr< betterSeedStruct> tmpDNAst = make_shared<betterSeedStruct>(tmp);

			bestDNA.push_back(tmpDNAst);

			oriKey.push_back(oriClKey);
			clusCnt = (int)bestDNA.size() - 1;
			seq2CI[segs2] = clusCnt;
		}
		else {//cluster exists
			cerr << "Found double id_ " << segs2 << endl << "Aborting.." << endl;
			exit(74);
		}
		addCnt++;
	}
	int siz = (int)bestDNA.size();
	//resize complete matrix to take these up
	
	for (uint i = 0; i < OTUmat.size(); i++){
		OTUmat[i].resize(siz, (matrixUnit)0);
	}
}

void UClinks::addSeedCounts(const vector<int>& clusters, shared_ptr<DNAunique> d) {
    for (int cluster : clusters) { add2OTUmat(d, cluster, (matrixUnit)clusters.size()); }
}

void UClinks::besterDNA(const vector<int>& curCLIDpre, shared_ptr<DNAunique> tdn1,
        Filters* fil, bool addCounts, uint64_t candidateAbundance) {
    if (curCLIDpre.empty()) { return; }
    if (addCounts) { addSeedCounts(curCLIDpre, tdn1); }
    // Multi-target (chimeric) assignments contribute counts, never seed candidates.
    if (curCLIDpre.size() != 1) { return; }
    const int curCLID = curCLIDpre[0];
    int TagIdx = -2;
	if (SeedsAreWritten){
		//delete dnaTemp1; if (dnaTemp2 != NULL){ delete dnaTemp2; }
		return;
	}
	//general routine, matching tdn found
	if (bestDNA[curCLID] == NULL) { //just fill with current DNA
		if (!seedSubclusters && fil->doReversePrimers() && !fil->checkYellowAndGreen(tdn1, CurSetPair, TagIdx, true, Filters::FilteringContext::Seeding)) {
			return;//delete dnaTemp1;
		}
		shared_ptr< betterSeedStruct> tmpDNAst = make_shared<betterSeedStruct>(tdn1);
		tmpDNAst->abundance = candidateAbundance;
		bestDNA[curCLID] = tmpDNAst;
		//fil->checkYellowAndGreen(tdn2, 1, TagIdx, true);
		//bestDNA2[curCLID] = tdn2;			
	}  	else { //already a candidate sequence? check who is better..
		bestDNA[curCLID]->whoIsBetter2(tdn1, candidateAbundance);
	}
}

bool betterSeedStruct::whoIsBetter2(shared_ptr<DNAunique> d1, uint64_t candidateAbundance){
	bool d1Merge = d1->getMerge() != nullptr;
	bool refMerge = tdn1->getMerge() != nullptr;
	if (d1Merge) {
		mergeCnt++; cummMergeLen += d1->getMerge()->length();
	} else { notMergeCnt++; }
	float mergeFreq = (float)mergeCnt / (mergeCnt + notMergeCnt);
	float avgMergeLen = mergeCnt > 0 ? float(cummMergeLen / mergeCnt) : 0.f;
	

	//fil->betterSeed(tdn1,  bestDNA[curCLID]->tdn1, bestDNA[curCLID]->bestPID, CurSetPair, checkBC) //bestLEN[curCLID], 
	bool formalWIB = whoIsBetter(d1, d1->getPair(), d1->getMerge(),
        tdn1, tdn1->getPair(), tdn1->getMerge(), bestPID, true, d1.get(),
        candidateAbundance ? candidateAbundance : d1->totalSum(),
        abundance ? abundance : tdn1->totalSum());

	bool replace(formalWIB);
	//had some other code here previously.. but formalWIB is now the main, generalized decision maker

	if (replace) {
		if (d1->getTempFloat()){
			bestPID = d1->getTempFloat();
		}
		abundance = candidateAbundance;
		len = d1->length(); 
		tdn1 = d1; 
		if (d1Merge) {
			mergLen = d1->getMerge()->length();
		}
	}

	return replace;
}

void UClinks::removeSampleID(string& w, const string &SEP) {
	size_t pos = w.find(SEP);
	if (pos != std::string::npos) {
		w = w.substr(pos + SEP.length());
	}
}
void UClinks::removeSampleID(string& w, const string &SEP, string & SMplID) {
	size_t pos = w.find(SEP);
	if (pos != std::string::npos) {
		SMplID = w.substr(0,pos);
		w = w.substr(pos + SEP.length());
	}
}
void UClinks::removeSizeStr(string& w) {
	size_t idx = w.find(";size=");
	if (idx != std::string::npos) {
		w = w.substr(0, idx);
	}
}
void UClinks::removeSeqStr(string& w) {
	size_t idx = w.find(";seqs=");
	if (idx != std::string::npos) {
		w = w.substr(0, idx);
	}
}
void UClinks::removeCentrStr(string& w) {
	size_t idx = w.find("centroid=");
	if (idx != std::string::npos) {
		w = w.substr(idx + 9);
	}
}

void UClinks::writeNewSeeds(OutputStreamer* MD, Filters* fil, 
			bool refSeeds, bool printLnk) {
	if (!RefDBmode && refSeeds){ return; }
	//ofstream O(outf.c_str());
	int paired = fil->isPaired();
	//MD->printStorage();
	ofstream links;
	uint st(0); uint to((uint)oriKey.size());
	//in case of refDB mode, need to start from point where refDBs were added in..
	if (refSeeds && RefDBmode && RefDBotuStart >= 0){
		cerr << "Writing ref DB sequences..";
		st = RefDBotuStart;
		//and also write a link between OTU_name and refSeq
		if (printLnk){
			string refLinkF = MD->leadOutFile() + ".lnks";
			links.open(refLinkF.c_str(), ios::out);
		}
	} else if (!refSeeds && RefDBmode && RefDBotuStart >= 0){
		cerr << "Writing new OTU seeds..";
		to = RefDBotuStart;
	}
	shared_ptr<DNAunique> d;


	
	for (uint i = st; i < to; i++) {
		//if (i == 6628){			int x = 0;		}
		//if check DNA was done before (remove rev primer), do it now
		if (bestDNA[i] == NULL) {
			//this should not happen on a regular lotus run
			MD->closeOutStreams();
			cerr << "No seed sequence found for DNA " << oriKey[i] << ". Aborting..\n"; // " << bestDNA[i]->getOldId()<<"(
			exit(54);
		}
		//fil->check(bestDNA[i],true);
		string newH = oriKey[i];
		d = bestDNA[i]->tdn1;
		shared_ptr<DerepRead> d2 = d->getPair();
		if (printLnk){
			string oriH = d->getShortId(); removeSizeStr(oriH);
			links << newH << "\t" << oriH << endl;
		}
		// Emit the selected consensus once, in the requested output format.
		if (paired == 2 && b_merge_pairs_optiSeed_ && MD->saveForWrite_merge(d, newH, 0, false)) {
			continue;
		}
		//takeOver id_ to allow for later linkup with cluster
		if (paired == 2 && !MD->singleSeedFasta()) {
			d->setPassed(true);
			d->setNewID(newH + ".1");


			int Cstream(100);
			if (d2 != nullptr) {
			    //sorted by pair1,2, quality yellow green, singleton (pair 1 2 3 4) etc
				d2->setPassed(true);
				d2->setNewID(newH + ".2");
				MD->saveForWrite(d, 1,-1, Cstream,false);
				MD->saveForWrite(d2, 2,-1, Cstream,false);
				MD->writeForWrite(d, 1, Cstream, d2, 2, Cstream);
				

			} else {
				MD->saveForWrite(d, 3,-1, Cstream,true);
			}
		} else {
			int Cstream(0);
			d->setNewID(newH);
			d->setPassed(true);
			MD->saveForWrite(d, 1,-1, Cstream,true);
		}
	}
	MD->closeOutStreams();
	SeedsAreWritten = true;
	if (printLnk){ links.close(); }
	cerr << "Done" << endl;
}
void UClinks::printStats(ostream& os){
	if (oriKey.size() == 0) {
		os << "No OTU's to re-seed\n";
		return;
	}
	uint numCls = 0;
	float avgQ(0.f), minQ(10000.f), maxQ(0.f),
		avgA(0.f), minA(10000.f), maxA(0.f),
		avgS(0.f), minS(10000.f), maxS(0.f);
	uint avgL = 0, minL(10000), maxL(0);
	uint merges(0);
	
	//uint div(0);// = (uint)oriKey.size();
	vector<uint> lengths;
	vector<float> quals,accums,sims;
	//oriKey.size() - not if there are only refs
	uint to((uint)oriKey.size());
	if (RefDBmode && RefDBotuStart >= 0){
		to = RefDBotuStart;
	}

	for (uint i = 0; i < to; i++){
		const shared_ptr<DNAunique>& candidate = bestDNA[i]->tdn1;
		if (candidate == nullptr) { continue; }
		shared_ptr<DerepRead> d = candidate;
		shared_ptr<DerepRead> d2(nullptr);
		if (candidate->getMerge() != nullptr) {
			d = candidate->getMerge();
			merges++;
		} else {
			d2 = candidate->getPair();
		}
		float curQ = d->getAvgQual();
		if (d2 != NULL) {
			curQ += d2->getAvgQual(); curQ /= 2.f;
		}
		uint curL = (uint) d->length();
		if (d2!= NULL) {
			curL += d2->length(); 
		}

		if (curL < minL){ minL = curL; }
		if (curL > maxL){ maxL = curL; }
		avgL += curL;
		lengths.push_back(curL);
		
		if (curQ < 1) {//no new Seed found, default seed
			continue;
		}
		float sc = bestDNA[i]->tdn1->getTempFloat(); // set on the candidate, not its merged product

		if (curQ < minQ){ minQ = curQ; }
		if (curQ > maxQ){ maxQ = curQ; }
		avgQ += curQ;

		float curA = (float)d->getAccumError();
		if (d2 != NULL) {
			curA += (float) d2->getAccumError(); 
		}

		quals.push_back(curQ);
		accums.push_back(curA); 
		if (curA < minA){ minA = curA; }
		if (curA > maxA){ maxA = curA; }
		avgA += curA;	

		sims.push_back(sc);
		if (sc < minS){ minS = sc; }
		if (sc > maxS){ maxS = sc; }
		avgS += sc;
		numCls++; 
	}
	os << "Found " << numCls << " seeds of " << oriKey.size() << " OTU's in " << uclines << " mappings." << endl;
	os << "Of these, " << merges << " were merged paired reads.\n";
    if (seedMergeAttempts_ > 0) {
        os << seedMergeSuccesses_ << "/" << seedMergeAttempts_ << " seed candidates had a paired-read overlap.\n";
    }

	//sort vectors
	std::sort(lengths.begin(), lengths.end());
	std::sort(quals.begin(), quals.end());
	std::sort(accums.begin(), accums.end());
	std::sort(sims.begin(), sims.end());
	os << "Stats of Seed sequences (0th/10th/50th/90th/100th) percentile:";
	if (lengths.size() > 0) { os << "\n     - Sequence Length :   " << calc_median2(lengths, 0.f) << "/" << calc_median2(lengths, 0.1f) << "/" << calc_median2(lengths, 0.5f) << "/" << calc_median2(lengths, 0.9f) << "/" << calc_median2(lengths, 1.f);; }
	if (quals.size() > 0) { os << "\n     - Quality :      " << calc_median2(quals, 0.f) << "/" << calc_median2(quals, 0.1f) << "/" << calc_median2(quals, 0.5f) << "/" << calc_median2(quals, 0.9f) << "/" << calc_median2(quals, 1.f);; }
	if (accums.size() > 0) { os << "\n     - Accum. Error : " << calc_median2(accums, 0.f) << "/" << calc_median2(accums, 0.1f) << "/" << calc_median2(accums, 0.5f) << "/" << calc_median2(accums, 0.9f) << "/" << calc_median2(accums, 1.f);;	}
	if (sims.size() > 0) {	os << "\n     - Sim2Consensus: " << calc_median2(sims, 0.f) << "/" << calc_median2(sims, 0.1f) << "/" << calc_median2(sims, 0.5f) << "/" << calc_median2(sims, 0.9f) << "/" << calc_median2(sims, 1.f);		}
	os << endl;
}
