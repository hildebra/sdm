// Direct sensitivity tests for the proposed seed rules.
#define SDM_99_CANDIDATES_NO_MAIN
#include "Test99Candidates.cpp"
bool intersect(const vector<EvalSeed>&a,const vector<EvalSeed>&b) {
    size_t i=0,j=0;
    while(i<a.size() && j<b.size()){if(a[i].key==b[j].key)return true;if(a[i].key<b[j].key)++i;else ++j;}
    return false;
}
bool signature_equal(const vector<EvalSeed>&a,const vector<EvalSeed>&b) {
    if(a.size()!=b.size())return false;
    for(size_t i=0;i<a.size();++i)if(a[i].key!=b[i].key || a[i].position!=b[i].position)return false;
    return true;
}
int main() {
    mt19937_64 rng(87367);
    cout<<"length,errors,trials,first4_all,first4_loose5,first4_budget,all_seed_budget,min50_any,min66_any,min50_identical_signature,mean_min50_seeds,mean_min66_seeds\n";
    for(size_t n:{100,150,200,250,300,400,1000}) {
        size_t e=n/100,trials=5000,strict=0,loose=0,soft=0,all=0,m50=0,m66=0,same=0,count50=0,count66=0;
        for(size_t t=0;t<trials;++t) {
            string a(n,'A');for(auto&c:a)c="ACGT"[rng()%4];string b=a;set<size_t> mutations;
            while(mutations.size()<e)mutations.insert(rng()%n);
            for(size_t p:mutations)b[p]="ACGT"[(base(b[p])+1+rng()%3)%4];
            auto grid=seeds(a,"grid15_scalar");size_t bad=0,bad4=0;bool looseOK=true;
            for(size_t j=0;j<grid.size();++j) {
                size_t p=grid[j].position;
                bool mismatch=a.compare(p,15,b,p,15)!=0;bad+=mismatch;
                if(j<4) {
                    bad4+=mismatch;
                    size_t start=p<5?0:p-5,end=min(p+5,n-15);bool found=false;
                    for(size_t x=start;x<=end;++x)if(a.compare(p,15,b,x,15)==0){found=true;break;}
                    looseOK &= found;
                }
            }
            strict+=bad4==0;loose+=looseOK;soft+=bad4<=e;all+=bad<=e;
            auto a50=seeds(a,"min15_w50"),b50=seeds(b,"min15_w50");
            auto a66=seeds(a,"min15_w66_heuristic"),b66=seeds(b,"min15_w66_heuristic");
            m50+=intersect(a50,b50);m66+=intersect(a66,b66);same+=signature_equal(a50,b50);
            count50+=a50.size();count66+=a66.size();
        }
        if(soft!=trials || all!=trials)throw runtime_error("Safe disjoint-seed rule rejected a true match");
        cout<<n<<','<<e<<','<<trials<<','<<strict<<','<<loose<<','<<soft<<','<<all<<','<<m50<<','<<m66<<','<<same<<','<<double(count50)/trials<<','<<double(count66)/trials<<endl;
    }
    bool adversarialFound=false, identicalSketchFound=false;
    // Produce an adversarial 200 bp / two-substitution minimizer miss.
    for(size_t t=0;t<100000;++t) {
        string a(200,'A');for(auto&c:a)c="ACGT"[rng()%4];auto selected=seeds(a,"min15_w66_heuristic");auto bypos=selected;
        sort(bypos.begin(),bypos.end(),[](auto&a,auto&b){return a.position<b.position;});
        vector<size_t> mutations;
        for(auto&seed:bypos)if(mutations.empty() || seed.position>mutations.back())mutations.push_back(seed.position+14);
        if(mutations.empty() || mutations.size()>2)continue;
        string b=a;for(size_t p:mutations)b[p]="ACGT"[(base(b[p])+1)%4];
        if(!intersect(selected,seeds(b,"min15_w66_heuristic"))) {
            adversarialFound=true;
            cerr<<"ADVERSARIAL w66 miss; mutations";for(auto p:mutations)cerr<<' '<<p;cerr<<"\nreference="<<a<<"\nquery="<<b<<'\n';break;
        }
    }
    // A whole-sketch hash cannot prove exact equality, even when positions are included.
    for(size_t t=0;t<100;++t) {
        string a(200,'A');for(auto&c:a)c="ACGT"[rng()%4];auto selected=seeds(a,"min15_w50");bool found=false;
        for(size_t p=0;p<a.size();++p) {
            string b=a;b[p]="ACGT"[(base(b[p])+1)%4];
            if(signature_equal(selected,seeds(b,"min15_w50"))) {
                cerr<<"IDENTICAL sketch, different sequences; mutation="<<p<<"\nreference="<<a<<"\nquery="<<b<<'\n';found=true;break;
            }
        }
        if(found){identicalSketchFound=true;break;}
    }
    if(!adversarialFound || !identicalSketchFound)throw runtime_error("Expected seed counterexample was not found");
    return 0;
}
