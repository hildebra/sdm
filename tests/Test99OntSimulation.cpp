// Rough full-length 16S-like Q20 simulation using the existing candidate benchmark.
// All reads are 1500 bp; errors are independent substitutions, without indels.
#define SDM_99_CANDIDATES_NO_MAIN
#include "Test99Candidates.cpp"

struct SimulatedRead {
    string sequence;
    size_t origin;
    size_t errors;
};

bool event(mt19937_64& rng, size_t per10000) {
    return rng()%10000 < per10000;
}
void mutate(char& c, mt19937_64& rng) {
    c="ACGT"[(base(c)+1+rng()%3)%4];
}
size_t mismatches(const string& a, const string& b) {
    size_t total=0;
    for(size_t p=0;p<a.size();++p)total+=a[p]!=b[p];
    return total;
}

// Explore streaming admission thresholds. Input fixture storage is outside timing;
// the candidate index only retains centroids. Reservoir counts are capacity estimates.
void compression_sweep(const vector<string>& reads,const vector<SimulatedRead>& simulated) {
    robin_hood::unordered_node_map<string,size_t> originBySequence;
    for(const auto& read:simulated)originBySequence[read.sequence]=read.origin;
    auto unique=items(reads);
    cout<<"repeat,identity_percent,seconds,reads,groups,candidates,index_entries,mixed_template_groups,minority_template_reads,mixed_lineage_groups,retain_5,retain_20,retain_50\n";
    for(size_t repeat=0;repeat<3;++repeat) {
        auto v=unique;mt19937_64 rng(2026091210ULL+repeat);
        shuffle(v.begin(),v.end(),rng);
        for(size_t budget:{10,15,20,25,30,40}) {
            auto start=Clock::now();auto result=cluster(v,"tile20_packed",budget);
            double seconds=chrono::duration<double>(Clock::now()-start).count();
            vector<map<size_t,size_t>> contents(v.size());
            for(size_t row=0;row<v.size();++row) {
                const size_t rep=result.assignment[row];
                if(mismatches(v[row].s,v[rep].s)>v[row].s.size()*budget/1000)
                    throw runtime_error("Assigned read exceeds streaming threshold");
                contents[rep][originBySequence.at(v[row].s)]+=v[row].count;
            }
            size_t mixed=0,minority=0,mixedLineages=0,retained5=0,retained20=0,retained50=0;
            for(const auto& histogram:contents) {
                if(histogram.empty())continue;
                size_t count=0,majority=0;set<size_t> lineages;
                for(const auto& [origin,n]:histogram){count+=n;majority=max(majority,n);lineages.insert(origin/10);}
                mixed+=histogram.size()>1;minority+=count-majority;mixedLineages+=lineages.size()>1;
                retained5+=min(size_t(5),count);retained20+=min(size_t(20),count);retained50+=min(size_t(50),count);
            }
            cout<<repeat<<','<<100.-double(budget)/10<<','<<seconds<<','<<reads.size()<<','<<result.groups<<','
                <<result.checks<<','<<result.entries<<','<<mixed<<','<<minority<<','<<mixedLineages<<','
                <<retained5<<','<<retained20<<','<<retained50<<endl;
        }
    }
}

int main(int argc,char** argv) {
    const size_t length=1500, lineageCount=5, variantsPerLineage=10;
    const size_t templateCount=lineageCount*variantsPerLineage;
    size_t readsPerTemplate=100;
    if(argc>=2)readsPerTemplate=stoul(argv[1]);
    if(argc>3 || (argc==3 && string(argv[2])!="--compression") || readsPerTemplate<2 || readsPerTemplate>10000) {
        cerr<<"Usage: "<<argv[0]<<" [reads-per-template: 2..10000; default 100] [--compression]\n";
        return 2;
    }
    mt19937_64 templateRng(2026091201ULL), readRng(2026091202ULL), diagnosticRng(2026091203ULL);
    string ancestral(length,'A');
    for(auto& c:ancestral)c="ACGT"[templateRng()%4];
    // Nine artificial 60 bp variable blocks, spread across a conserved backbone.
    // This mimics conserved/variable structure, not a biological 16S sequence model.
    vector<bool> variable(length,false);
    for(size_t block=0;block<9;++block)
        for(size_t p=80+block*155;p<140+block*155;++p)variable[p]=true;
    vector<string> references;
    for(size_t lineage=0;lineage<lineageCount;++lineage) {
        string parent=ancestral;
        for(size_t p=0;p<length;++p)if(variable[p] && event(templateRng,2000))mutate(parent[p],templateRng);
        for(size_t variant=0;variant<variantsPerLineage;++variant) {
            string reference=parent;
            for(size_t p=0;p<length;++p)if(variable[p] && event(templateRng,200))mutate(reference[p],templateRng);
            references.push_back(move(reference));
        }
    }
    vector<SimulatedRead> simulated;
    vector<string> reads;
    uint64_t totalErrors=0;
    size_t withinReferenceBudget=0,minErrors=length,maxErrors=0;
    for(size_t origin=0;origin<references.size();++origin) {
        for(size_t copy=0;copy<readsPerTemplate;++copy) {
            string sequence=references[origin];size_t errors=0;
            for(auto& c:sequence)if(event(readRng,100)){mutate(c,readRng);++errors;}
            if(mismatches(sequence,references[origin])!=errors)throw runtime_error("Error simulation mismatch");
            totalErrors+=errors;withinReferenceBudget+=errors<=length/100;
            minErrors=min(minErrors,errors);maxErrors=max(maxErrors,errors);
            reads.push_back(sequence);simulated.push_back({move(sequence),origin,errors});
        }
    }
    double withinDivergence=0,betweenDivergence=0;size_t withinPairs=0,betweenPairs=0;
    for(size_t a=0;a<references.size();++a)for(size_t b=a+1;b<references.size();++b) {
        size_t errors=mismatches(references[a],references[b]);
        if(a/variantsPerLineage==b/variantsPerLineage){withinDivergence+=errors;++withinPairs;}
        else {betweenDivergence+=errors;++betweenPairs;}
    }
    size_t pairTrials=10000,pairMatches=0,pairFirst4=0;
    uint64_t pairErrors=0;
    ReadMerger verifier;
    for(size_t trial=0;trial<pairTrials;++trial) {
        size_t origin=diagnosticRng()%templateCount;
        size_t a=diagnosticRng()%readsPerTemplate,b=diagnosticRng()%(readsPerTemplate-1);
        if(b>=a)++b;
        const string& qa=simulated[origin*readsPerTemplate+a].sequence;
        const string& qb=simulated[origin*readsPerTemplate+b].sequence;
        size_t errors=mismatches(qa,qb);pairErrors+=errors;
        bool accepted=errors<=length/100;
        Item q{qa,1,{}},r{qb,1,{}};pack(q);pack(r);
        if(packed_match(q,r,length/100)!=accepted ||
           (verifier.percentIdentity(qa.data(),qb.data(),int(length),int(length/100))>=0)!=accepted)
            throw runtime_error("1500 bp verifier disagrees with direct mismatch count");
        if(accepted) {
            ++pairMatches;bool first4=true;
            for(size_t p:{0,33,66,100})first4 &= qa.compare(p,15,qb,p,15)==0;
            pairFirst4+=first4;
        }
    }
    cerr<<setprecision(10)
        <<"simulation: seed_templates=2026091201 seed_reads=2026091202 length="<<length
        <<" templates="<<templateCount<<" lineages="<<lineageCount<<" reads="<<reads.size()<<'\n'
        <<"errors: total="<<totalErrors<<" mean="<<double(totalErrors)/reads.size()
        <<" rate="<<double(totalErrors)/(length*reads.size())<<" min="<<minErrors<<" max="<<maxErrors
        <<" within_99pct_of_true_template="<<withinReferenceBudget<<'/'<<reads.size()<<'\n'
        <<"template_divergence: within_lineage="<<withinDivergence/(withinPairs*length)
        <<" between_lineages="<<betweenDivergence/(betweenPairs*length)<<'\n'
        <<"same_template_read_pairs: trials="<<pairTrials<<" mean_mismatches="<<double(pairErrors)/pairTrials
        <<" within_99pct="<<pairMatches<<" first4_also_pass="<<pairFirst4<<'\n';

    if(argc==3) {compression_sweep(reads,simulated);return 0;}

    auto v=items(reads);
    auto start=Clock::now();auto truth=cluster(v,"brute");
    cerr<<"exhaustive_reference: seconds="<<chrono::duration<double>(Clock::now()-start).count()
        <<" unique="<<v.size()<<" groups="<<truth.groups<<" candidates="<<truth.checks<<'\n';
    auto indexed=cluster(v,"tile20_packed");
    if(indexed.assignment!=truth.assignment)throw runtime_error("Full simulation differs from exhaustive clustering");
    robin_hood::unordered_node_map<string,size_t> originBySequence;
    for(auto& read:simulated)originBySequence[read.sequence]=read.origin;
    size_t crossTemplateAssignments=0;
    vector<set<size_t>> templateClusters(templateCount);
    for(size_t i=0;i<v.size();++i) {
        size_t origin=originBySequence.at(v[i].s),representative=indexed.assignment[i];
        crossTemplateAssignments+=origin!=originBySequence.at(v[representative].s);
        templateClusters[origin].insert(representative);
    }
    size_t minClusters=readsPerTemplate,maxClusters=0;
    for(auto& groups:templateClusters){minClusters=min(minClusters,groups.size());maxClusters=max(maxClusters,groups.size());}
    cerr<<"cluster_diagnostics: cross_template_assignments="<<crossTemplateAssignments
        <<" clusters_per_template_min="<<minClusters<<" max="<<maxClusters<<'\n';

    const vector<string> methods={"tile20_scalar","tile20_packed","grid15_scalar","grid15_packed",
        "grid15_vote","grid15_first4_strict","min15_w50_packed",
        "min15_w50_position_packed","min15_w66_position_packed"};
    cout<<"dataset,repeat,method,seconds,unique,groups,candidates,seed_rejections,index_entries,query_seeds,assignment_differences\n";
    run("ont_q20_16s_1500",reads,methods);
    return 0;
}
