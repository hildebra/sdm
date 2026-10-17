// Standalone candidate-index experiment; not the production clustering path.
#include "ReadMerger.h"
#include "include/robin_hood.h"
#include <bit>
#include <chrono>
#include <deque>
#include <iomanip>
#include <numeric>
#include <optional>
#include <random>
#include <set>
using namespace std;
using Clock = chrono::steady_clock;
struct Item { string s; uint64_t count; vector<uint64_t> packed; bool canonical=true; };
struct EvalSeed { size_t position; uint64_t key; };
struct Result { vector<size_t> assignment; size_t groups=0, checks=0, rejected=0, entries=0, seeds=0; };
int base(char c) { switch(c) {case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:return -1;} }
void pack(Item& a) {
    a.packed.assign((a.s.size()+31)/32,0);
    for(size_t i=0;i<a.s.size();++i) {
        int b=base(a.s[i]); if(b<0) {a.canonical=false;b=0;}
        a.packed[i/32] |= uint64_t(b)<<((i%32)*2);
    }
}
optional<uint64_t> word(const string& s,size_t pos,size_t k) {
    uint64_t out=0;
    for(size_t j=0;j<k;++j) {int b=base(s[pos+j]);if(b<0)return {};out=(out<<2)|b;}
    return out;
}
uint64_t mix(uint64_t x) {
    x ^= x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);
}
vector<EvalSeed> seeds(const string& s,const string& method) {
    vector<EvalSeed> out;
    if(method.starts_with("min15_")) {
        if(s.size()<15)return out;
        const size_t n=s.size()-14,w=min(size_t(method.starts_with("min15_w50")?50:66),n);
        vector<optional<uint64_t>> words(n);vector<uint64_t> hashes(n);
        uint64_t rolling=0;size_t valid=0;
        for(size_t i=0;i<s.size();++i) {
            int b=base(s[i]);valid=b<0?0:valid+1;
            rolling=((rolling<<2)|uint64_t(b<0?0:b))&((uint64_t(1)<<30)-1);
            if(i>=14){size_t pos=i-14;if(valid>=15)words[pos]=rolling;hashes[pos]=words[pos]?mix(*words[pos]):~uint64_t(0);}
        }
        deque<size_t> q;
        size_t previous=~size_t(0);
        for(size_t i=0;i<n;++i) {
            while(!q.empty() && q.front()+w<=i)q.pop_front();
            while(!q.empty() && hashes[q.back()]>hashes[i])q.pop_back();
            q.push_back(i);
            if(i+1>=w && q.front()!=previous && words[q.front()]) {
                previous=q.front();
                const uint64_t key = method.find("position") == string::npos ? *words[previous] :
                    (uint64_t(previous)<<30)|*words[previous];
                out.push_back({previous,key});
            }
        }
        // One posting per representative/token, regardless of repeat count.
        sort(out.begin(),out.end(),[](auto&a,auto&b){return a.key<b.key;});
        out.erase(unique(out.begin(),out.end(),[](auto&a,auto&b){return a.key==b.key;}),out.end());
    } else {
        bool tile20=method.starts_with("tile20");size_t k=tile20?20:15;
        for(size_t slot=0;;++slot) {
            size_t pos=tile20?slot*20:(slot/3)*100+(slot%3)*33;
            if(pos+k>s.size())break;
            auto value=word(s,pos,k);
            if(value)out.push_back({pos,(uint64_t(slot)<<40)|*value});
        }
    }
    return out;
}
bool packed_match(const Item& q,const Item& r,size_t e) {
    size_t mismatches=0,n=q.s.size();
    for(size_t j=0;j<q.packed.size();++j) {
        uint64_t x=q.packed[j]^r.packed[j];
        if(j+1==q.packed.size() && n%32)x&=(uint64_t(1)<<(2*(n%32)))-1;
        mismatches+=popcount((x|(x>>1))&0x5555555555555555ULL);
        if(mismatches>e)return false;
    }
    return true;
}
vector<Item> items(const vector<string>& reads) {
    robin_hood::unordered_node_map<string,uint64_t> count;
    for(auto&s:reads)++count[s];
    vector<Item> out;out.reserve(count.size());
    for(auto&[s,n]:count){out.push_back({s,n,{}});pack(out.back());}
    sort(out.begin(),out.end(),[](auto&a,auto&b){
        if(a.s.size()!=b.s.size())return a.s.size()>b.s.size();
        if(a.count!=b.count)return a.count>b.count;
        return a.s<b.s;
    });
    return out;
}
Result cluster(const vector<Item>& v,const string& method,size_t errorsPerThousand=10) {
    if(errorsPerThousand>1000)throw invalid_argument("Invalid mismatch rate");
    using Index=robin_hood::unordered_node_map<uint64_t,vector<size_t>>;
    Index index;index.reserve(v.size()*4);
    vector<size_t> reps;
    Result answer;ReadMerger verifier;
    bool packed=method.ends_with("packed"),strict=method=="grid15_first4_strict",vote=method=="grid15_vote";
    for(size_t row=0;row<v.size();++row) {
        const Item& query=v[row];size_t e=query.s.size()*errorsPerThousand/1000;
        auto querySeeds=seeds(query.s,method);answer.seeds+=querySeeds.size();
        bool mini=method.starts_with("min15_");
        bool brute=method=="brute" || (mini ? querySeeds.empty() : querySeeds.size()<e+1);
        // w=50 guarantees a common minimizer if L >= (e+1)*(k+w-1).
        // Use exhaustive lookup for shorter/ambiguous queries in this variant.
        if(method.starts_with("min15_w50") && (query.s.size()<64*(e+1) || !query.canonical))brute=true;
        vector<const vector<size_t>*> lists;
        if(!brute) {
            for(auto&seed:querySeeds) {
                auto found=index.find(seed.key);
                lists.push_back(found==index.end()?nullptr:&found->second);
            }
            sort(lists.begin(),lists.end(),[](auto a,auto b){return (a?a->size():0)<(b?b->size():0);});
            if(!mini)lists.resize(e+1);
        }
        vector<size_t> cursor(lists.size(),0);
        size_t winner=reps.size(),bruteCursor=0;
        while(true) {
            size_t next=reps.size();
            if(brute) {if(bruteCursor<reps.size())next=bruteCursor++;}
            else for(size_t j=0;j<lists.size();++j)if(lists[j] && cursor[j]<lists[j]->size())next=min(next,(*lists[j])[cursor[j]]);
            if(next==reps.size())break;
            ++answer.checks;
            const Item& candidate=v[reps[next]];
            bool seedPass=true;
            if(strict || vote) {
                size_t bad=0,limit=strict?min(size_t(4),querySeeds.size()):querySeeds.size();
                for(size_t j=0;j<limit;++j) {
                    if(memcmp(query.s.data()+querySeeds[j].position,candidate.s.data()+querySeeds[j].position,15))++bad;
                    if(bad>(strict?0:e)){seedPass=false;break;}
                }
            }
            if(!seedPass)++answer.rejected;
            bool match=seedPass && (packed && query.canonical && candidate.canonical ? packed_match(query,candidate,e) :
                verifier.percentIdentity(query.s.data(),candidate.s.data(),int(query.s.size()),int(e))>=0.);
            if(match){winner=next;break;}
            if(!brute)for(size_t j=0;j<lists.size();++j)if(lists[j] && cursor[j]<lists[j]->size() && (*lists[j])[cursor[j]]==next)++cursor[j];
        }
        if(winner==reps.size()) {
            reps.push_back(row);
            if(method!="brute")for(auto&seed:querySeeds){index[seed.key].push_back(winner);++answer.entries;}
        }
        answer.assignment.push_back(reps[winner]);
    }
    answer.groups=reps.size();return answer;
}
void validate() {
    mt19937_64 rng(9032);vector<string> reads;
    for(size_t family=0;family<24;++family) {
        string reference(620,'A');for(auto&c:reference)c="ACGT"[rng()%4];
        for(size_t n:{14,15,31,64,99,100,150,199,200,250,299,300,399,400,500}) {
            reads.push_back(reference.substr(0,n));
            for(size_t errors=1;errors<=n/100+1;++errors) {
                string s=reference.substr(0,n);
                for(size_t j=0;j<errors;++j){size_t p=(j*33+7)%n;s[p]="ACGT"[(base(s[p])+1)%4];}
                reads.push_back(s);
            }
        }
    }
    reads.push_back("ACGTNNNNACGTACGT");reads.push_back("ACGTNNNNACGTACGTA");
    auto v=items(reads);auto truth=cluster(v,"brute");
    for(auto method:{"tile20_scalar","tile20_packed","grid15_scalar","grid15_packed","grid15_vote","min15_w50","min15_w50_position","min15_w50_position_packed"}) {
        auto indexed=cluster(v,method);
        if(indexed.assignment!=truth.assignment)throw runtime_error(string(method)+" differs from brute-force clustering");
    }
    // Differential verifier check, including tail masks and mutation boundaries.
    ReadMerger verifier;
    for(size_t i=0;i<5000;++i) {
        size_t n=1+rng()%600;string s(n+90,'A');for(auto&c:s)c="ACGT"[rng()%4];
        Item r{s,1,{}},q{s.substr(0,n),1,{}};
        for(size_t j=0;j<i%12;++j){size_t p=rng()%n;q.s[p]="ACGT"[rng()%4];}
        pack(r);pack(q);
        bool expected=verifier.percentIdentity(q.s.data(),r.s.data(),int(n),int(n/100))>=0;
        if(packed_match(q,r,n/100)!=expected)throw runtime_error("packed mismatch count disagrees");
    }
    cerr<<"PASS "<<v.size()<<" boundary/ambiguity fixtures against exhaustive clustering; 5000 packed-verifier comparisons\n";
}
void run(const string& name,const vector<string>& reads,const vector<string>& methods) {
    auto start=Clock::now();auto v=items(reads);
    cerr<<name<<": "<<reads.size()<<" reads, "<<v.size()<<" unique; exact collapse + packing "<<chrono::duration<double>(Clock::now()-start).count()<<" s\n";
    auto expected=cluster(v,"tile20_scalar").assignment;
    for(size_t repeat=0;repeat<4;++repeat)for(auto&method:methods) {
        start=Clock::now();auto r=cluster(v,method);double seconds=chrono::duration<double>(Clock::now()-start).count();
        size_t differences=0;for(size_t i=0;i<v.size();++i)differences+=r.assignment[i]!=expected[i];
        if(method!="grid15_first4_strict" && method!="min15_w66_heuristic" && differences)throw runtime_error("lossless index assignment mismatch");
        cout<<name<<','<<repeat<<','<<method<<','<<seconds<<','<<v.size()<<','<<r.groups<<','<<r.checks<<','<<r.rejected<<','<<r.entries<<','<<r.seeds<<','<<differences<<endl;
    }
}
#ifndef SDM_99_CANDIDATES_NO_MAIN
int main(int argc,char**argv) {
    if(argc==2 && string(argv[1])=="--self-test") {validate();return 0;}
    if(argc<2 || argc>3) {
        cerr<<"Usage: "<<argv[0]<<" keys.tsv [comma-separated methods]\n       "<<argv[0]<<" --self-test\n";
        return 2;
    }
    validate();ifstream input(argv[1]);vector<string> reads;string line;
    vector<string> methods={"tile20_scalar","tile20_packed","grid15_scalar","grid15_packed","grid15_vote","grid15_first4_strict","min15_w50","min15_w66_heuristic"};
    if(argc==3) {
        const vector<string> allowed={"tile20_scalar","tile20_packed","grid15_scalar","grid15_packed","grid15_vote","grid15_first4_strict","min15_w50","min15_w50_position","min15_w50_position_packed","min15_w66_heuristic"};
        methods.clear();istringstream names(argv[2]);string method;
        while(getline(names,method,',')) {
            if(find(allowed.begin(),allowed.end(),method)==allowed.end())throw runtime_error("Unknown method: "+method);
            methods.push_back(method);
        }
        if(methods.empty())throw runtime_error("No benchmark methods selected");
    }
    while(getline(input,line))if(line.starts_with("1\t"))reads.push_back(line.substr(2));
    if(reads.empty())throw runtime_error("No captured passing keys");
    cout<<"dataset,repeat,method,seconds,unique,groups,candidates,seed_rejections,index_entries,query_seeds,assignment_differences\n";
    run("local_full",reads,methods);
    for(auto&s:reads)if(s.size()>200)s.resize(200);
    run("local_cap200",reads,methods);
    mt19937_64 rng(555);reads.clear();string common(115,'A');for(auto&c:common)c="ACGT"[rng()%4];
    for(size_t i=0;i<5000;++i){string s=common;s.resize(300,'A');for(size_t p=115;p<s.size();++p)s[p]="ACGT"[rng()%4];reads.push_back(s);}
    run("conserved_first115",reads,methods);
    return 0;
}

#endif
