//#define DBGMOD

#ifdef DBGMOD
//#define SHOWBLK
//#define SHOWSOL
//#define SHOWRECUR
#endif

#include <fstream>
#include "EnuBundle.h"
#include "Utility.h"
using namespace std;

/**
TODO: reverse[] is not supported.
*/
int EnuBundle::readBinaryGraph(const char* filepath) {
	FILE *f = Utility::open_file(filepath, "rb");
	ui tt;
	fread(&tt, sizeof(ui), 1, f);
	if (tt != sizeof(ui)) {
		printf("sizeof unsigned int is different: file %u, machine %lu\n", tt, sizeof(ui));
	}
	fread(&n, sizeof(ui), 1, f);
	fread(&m, sizeof(ui), 1, f);
	printf("n=%u, m=%u (undirected)\n", n, m);
	ui* degree = new ui[n];
	fread(degree, sizeof(ui), n, f);

	if (pstart != nullptr) delete[] pstart;
	pstart = new ui[n + 1];
	if (edges != nullptr) delete[] edges;
	edges = new ui[m];

	//if (reverse != nullptr) delete[] reverse;
	//reverse = new ui[m];

	pstart[0] = 0;
	for (ui i = 0; i < n; i++) {
		if (degree[i] > 0) fread(edges + pstart[i], sizeof(ui), degree[i], f);
		pstart[i + 1] = pstart[i] + degree[i];
	}
	fclose(f);
	delete[] degree;
	
	return 0;
}

int EnuBundle::degeneracyOrder(ui *seq, ui *core, ui* pos) {
	ui *id_s = seq, *degree = core;
	for (ui i = 0; i < n; i++) {
		id_s[i] = i;
		degree[i] = pstart[i + 1] - pstart[i];
	}

	ui max_core = 0;
	ListLinearHeap *linear_heap = new ListLinearHeap(n, n - 1);
	linear_heap->init(n, n - 1, id_s, degree);
	memset(core, 0, sizeof(ui)*n);
	for (ui i = 0; i < n; i++) {
		ui u, key;
		linear_heap->pop_min(u, key);
		if (key > max_core)
			max_core = key;
		seq[i] = u;
		core[u] = max_core;
		pos[u] = i;

		for (ui j = pstart[u]; j < pstart[u + 1]; j++)
			if (core[edges[j]] == 0) {
				linear_heap->decrement(edges[j]);
			}
	}
	delete linear_heap;

	return max_core;
}

ui EnuBundle::markBlock1(ui v, ui* adress) {
	memset(dist, (ui)0, sizeof(ui)* n);
	dist[v] = 1;
	ui szblk = 1;
	//vector<ui> nei1s;
	szc1 = 0;
	//vector<ui> blkset;
	szc2 = 0; //all vertex in the block
	for (ui j = pstart[v]; j < pstart[v + 1]; j++) {
		ui nei = edges[j];
		if (!dist[nei] && dpos[nei] > dpos[v]) {
			dist[nei] = 1;	//N(v)
			//blkset.emplace_back(nei);
			//nei1s.emplace_back(nei);
			cache1[szc1++] = nei;
			cache2[szc2++] = nei;
			szblk++;
		}
	}
	for (ui i = 0; i < szc1; i++) {
		ui u = cache1[i];
		//for (auto u:nei1s){
		for (ui k = pstart[u]; k < pstart[u + 1]; k++) {
			ui nei2 = edges[k];
			if (!dist[nei2] && dpos[nei2] > dpos[v]) {
				dist[nei2] = 2; //N2(v)
				//blkset.emplace_back(nei2);
				cache2[szc2++] = nei2;
				szblk++;
			}
		}
	}
	/*size pruning*/
	for (ui i = 0; i < szc2; i++) {
		ui u = cache2[i];
			//set<ui> vnei(edges + pstart[v], edges + pstart[v + 1]);
			//caclculate intersections
		ui cnt = 0;
		ui p1 = pstart[v], p2 = pstart[u];
		while (p1 != pstart[v + 1] && p2 != pstart[u + 1]) {
			if (edges[p1] < edges[p2]) {
				p1++;
			}
			else if (edges[p1] > edges[p2]) {
				p2++;
			}
			else {//equal
				if (dist[edges[p1]] ) {
					cnt++;
				}
				p1++, p2++;
			}
		}
		if (dist[u] == 1) {
			if (cnt + 2 * k < lb) {
				dist[u] = 0;
				//canprune = 1;
				szblk--;
			}
		}
		else {
			if (cnt + 2 * k - 2 < lb) {
				dist[u] = 0;
				//canprune = 1;
				szblk--;
			}
		}

	}
	ui cnt = 0;
	adress[cnt++] = v;
	for (ui i = 0; i < szc2; i++) {
		if (cache2[i]!=v &&dist[cache2[i]] != 0 )
			adress[cnt++] = cache2[i];
	}
	// assert(cnt == szblk);
	return szblk;
}


ui EnuBundle::markBlock2(ui v, ui* adress) {

	memset(dist, (ui)n, sizeof(ui)*n);
	memset(common, (ui)0, sizeof(ui)*n);
	
	szc1 = 0;
	for (ui i = pstart[v]; i < pstart[v + 1]; i++) {
		ui u = edges[i];		
		if (dpos[u] > dpos[v]) {
			dist[u] = 1;
			cache1[szc1++] = u;
		}
	}

	szc2 = 0;
	for (ui i = 0; i < szc1; i++) {
		ui u = cache1[i];
		for (ui j = pstart[u]; j < pstart[u + 1]; j++) {
			ui w = edges[j];
			if (w != v && dpos[w] > dpos[v]) {
				if (dist[w] == 1) {
					common[u]++;
				}
				else if (dist[w] == 2) {
					common[w]++;
				}
				else{
					common[w]++;					
					dist[w] = 2;
					cache2[szc2++] = w;
				}
			}
		}
	}
	ui szblk = 0;
	adress[szblk++] = v;
	for (ui i = 0; i < szc1; i++) {
		if (common[cache1[i]] + 2 * k >= lb) 
			adress[szblk++] = cache1[i];
	}
	for (ui i = 0; i < szc2; i++) {
		if (common[cache2[i]] + 2 * k - 2 >= lb)
			adress[szblk++] = cache2[i];
	}

	return szblk;
}

inline ui EnuBundle::isInBlock(ui vtx) {
	return nID[vtx] < bn;
}
int EnuBundle::buildBlock(ui v, ui *blk, ui sz) {
	bvtx = v; // the start vertex
	bn = sz;	// the size

	memset(bID, n, sizeof(ui)*n);
	memset(nID, n, sizeof(ui)*n);	
	for (ui i = 0; i < sz; i++) {
		bID[i] = blk[i];
		nID[blk[i]] = i;
	}

	//initialize the subgraph	
	/*if (bstart != nullptr) delete[] bstart;
	if (bedges != nullptr) delete[] bedges;	
	bstart = new ui[bn+1];
	bedges = new ui[bn*bn];
	*/

	bm = 0;
	for (ui j = 0; j < bn; j++) {
		int oruid = bID[j];
		badc[j].reinit(bn);
		//binv[j].init(bn);
		bstart[j] = bm;
		for (ui k = pstart[oruid]; k < pstart[oruid + 1]; k++) {
			if (isInBlock(edges[k])) {
				bedges[bm++] = nID[edges[k]];
				badc[j].set(nID[edges[k]]);
			}
		}
		//binv[j].flip();
		//sort(bedges + bstart[j], bedges + bm); //we can optimize here by sorting the vertex
	}
	bstart[bn] = bm;
	return 1;

}

inline void EnuBundle::CandToP(const ui &u) {
	// assert(Cand.contains(u));
	P.add(u);
	Cand.remove(u);

	//update neiInP
	for (ui i =  bstart[u]; i < bstart[u + 1]; i++) {
		ui &nei = bedges[i];
		++neiInP[nei];
	}	

}

inline void EnuBundle::PToCand(const ui &u) {
	// assert(P.contains(u));
	P.remove(u);
	Cand.add(u);

	//update neiInP	
	for (ui i = bstart[u]; i < bstart[u + 1]; i++) {
		ui &nei = bedges[i];
		--neiInP[nei];		
	}
}

void EnuBundle::removeFrCand(ui u) {
	// assert(Cand.contains(u));
	Cand.remove(u);
	for (ui i = bstart[u]; i < bstart[u + 1]; i++) {
		ui nei = bedges[i];
		neiInG[nei]--;
	}
}

void EnuBundle::addToCand(ui u) {
	// assert(!Cand.contains(u)); 
	Cand.add(u);
	for (ui i = bstart[u]; i < bstart[u + 1]; i++) {
		ui nei = bedges[i];
		neiInG[nei]++;
	}
}

int EnuBundle::canMoveToP(ui u) {
	//// assert(Cand.contains(u));
	if (neiInP[u] + k < P.getSize() + 1) {
		return 0;
	}
	for (ui i = 0; i < P.getSize(); i++) {
		ui v = P.get(i);
		// assert(neiInP[v] + k >= P.getSize());
		if (neiInP[v] + k == P.getSize() && (!badc[v].test(u))) {
			return 0;
		}
	}
	return 1;	
}
/*
void EnuBundle::updateSatSet() {
	satInP.clear();
	for (ui i = 0; i < P.getSize(); i++) {
		ui u = P.get(i);
		if (neiInP[u] + k == P.getSize())
			satInP.add(u);
	}
}*/

/**
lst1 and lst2 must be sorted
*/
ui EnuBundle::interSection(ui *lst1, ui sz1, ui *lst2, ui sz2, ui* dest) {
	ui i = 0, j = 0;
	ui szdest = 0;
	while (i < sz1 && j < sz2) {
		if (lst1[i] < lst2[j]) ++i;
		else if (lst1[i] > lst2[j])++j;
		else {
			dest[szdest++] = lst1[i];
			i++, j++;
		}
	}
	return szdest;
}

int EnuBundle::isGlobalMaximal() {
	//vector<ui> candidates;
	ui init = 0;
	szc1= 0; //saturated vertices
	//updateSatSet();
	//if (!satInP.empty()) {
	for (ui i = 0; i < P.getSize(); i++) {
		if (neiInP[P.get(i)] + k == P.getSize()) {
			cache1[szc1++] = P.get(i);
		}
	}
	if (szc1 != 0) {		
		ui orsat0 = bID[cache1[0]];
		szc2 = pstart[orsat0 + 1] - pstart[orsat0];
		copy(edges + pstart[orsat0], edges + pstart[orsat0 + 1], cache2);		
		for (ui i = 1; i < szc1; i++) {
			ui orsat_i = bID[cache1[i]];
			szc2 = interSection(edges + pstart[orsat_i], pstart[orsat_i + 1] - pstart[orsat_i],
				cache2, szc2, cache2);
		}
		for (ui i = 0; i < szc2; i++) {
			ui oru = cache2[i];
			if (dpos[oru] < dpos[bvtx] && canGloablaAdd(oru)) {
				return 0;
			}
		}
		return 1;
	}
	else {
		for (ui i = 0; i < dpos[bvtx]; i++) {
			ui u = dseq[i];
			if (core[u] + k >= P.getSize() + 1 &&  canGloablaAdd(u))
				return 0;
		}
		return 1;
	}
}
//Prerequist, oru is adjacent to all saturated vertices in P.
//Check the  |P\cup N(oru)| > |P|+1-k
int EnuBundle::canGloablaAdd(ui oru) {
	//int scon = 0;
	int pcon = 0;	
	for (ui j = pstart[oru]; j < pstart[oru + 1]; j++) {
		ui nei = edges[j];
		if (isInBlock(nei) && P.contains(nID[nei])) {
				pcon++;
		}
	}
	if (pcon + k >= P.getSize() + 1)
		return 1; //can add
	return 0;	//can not add
}

/**
A stronger branch rules
Move at most szmax vertices from doing to P.
*/
void EnuBundle::multiRecurSearch(vector<ui> &doing, ui szmax) {
	vector<ui> rcand;
	vector<ui> rexcl;
	// assert(szmax < doing.size());
	// assert(szmax > 0 && szmax < k);


	if (Cand.getSize() + P.getSize() < lb) { // prune
		return;
	}
#ifdef SHOWRECUR
	printf("At most %u from Doing: ", szmax);
	for (ui u : doing) {
		printf("%u ", u);
	}
	printf("\n");
#endif
	ui idx = 0;
	ui stop = 0;
	while(!stop && idx <= szmax){ // szmax +1 branches
		if (idx > 0) {
			ui ubr = doing[idx - 1];			
			CandToP(ubr);
#ifdef SHOWRECUR
			printf("Add %u\n", ubr);
#endif
			//vector<ui> tmpcand;
			szc1 = 0;
			//update candidate,
			//WARNING: iterator of Cand(Excl) will be invalid if it is deleted
			for (ui j = 0; j < Cand.getSize(); j++) {
				ui u = Cand.get(j);
				if (!canMoveToP(u)) {
					//tmpcand.push_back(u);					
					cache1[szc1++] = u;
				}
			}
			for (ui i = 0; i < szc1; i++) {
				removeFrCand(cache1[i]);
			}
			rcand.insert(rcand.end(), cache1, cache1+szc1);
			//update excl
			szc2 = 0;
			for (ui j = 0; j < Excl.getSize(); j++) {
				ui u = Excl.get(j);
				if (!canMoveToP(u)) {
					//tmpexcl.push_back(u);
					cache2[szc2++] = u;
				}
			}
			for (ui i = 0; i < szc2; i++) {
				Excl.remove(cache2[i]);
			}
			rexcl.insert(rexcl.end(), cache2, cache2+szc2);
		}
		
		//remove the last vertex all the remaining vertex;
		if (idx < szmax) {// the 0 to sz-1 branches
			if (!Cand.contains(doing[idx])) {// 
				//doing[idx] is reduced when adding doing[0->idx-1]. indicating
				//backtrack idx since we only add 0-(idx-1) vertices
				Excl.add(doing[idx]);
#ifdef SHOWRECUR
				printf("Remove %u, alread reduced\n", doing[idx]);
#endif
				branch();
				Excl.remove(doing[idx]);
				stop = 1;
			}
			else {
				removeFrCand(doing[idx]);
				Excl.add(doing[idx]);
#ifdef SHOWRECUR
				printf("Remove %u \n", doing[idx]);
#endif
				branch();
				addToCand(doing[idx]);
				Excl.remove(doing[idx]);
			}
		}
		else { // idx== szmax, the last branch
			vector<ui> doreduced;
			for (ui j = idx; j < doing.size(); j++) {
				if (Cand.contains(doing[j])) { //not reduced, we need to recover
					removeFrCand(doing[j]);
					doreduced.emplace_back(doing[j]);
#ifdef SHOWRECUR
					printf("Remove %u \n", doing[j]);
#endif
				}
				else {
#ifdef SHOWRECUR
					printf("Remove %u, alread reduced\n", doing[j]);
#endif
				}
				//Excl.add(doing[j]);
			}
			branch();

			//recover
			//for (ui j = idx; j < doing.size(); j++) {
			//	Excl.remove(doing[j]);
			//}
			for (auto u : doreduced)
				addToCand(u);			
		}
		idx++;
	}
	for (ui i = 0; i < idx-1; i++) {
		PToCand(doing[i]);
	}
	//recover
	for (auto u : rcand) {
		addToCand(u);
	}
	for (auto u : rexcl) {
		Excl.add(u);
	}
}

/**
TODO:
move start into solution and continue the search
*/
void EnuBundle::multiRecurSearch(ui start) {

}
/*
* Take a vertex start from cand and continue the search
*/
void EnuBundle::recurSearch(ui start) {
	//reduce
	// assert(Cand.contains(start));
	//// assert(canMoveToP(start));
	vector<ui> rcand;
	vector<ui> rexcl;

	if (Cand.getSize() + P.getSize() < lb) {
		return;
	}	
	CandToP(start);
	//update candidate
	for (ui i = 0; i < Cand.getSize(); i++) {
		ui u = Cand.get(i);
		if (!canMoveToP(u)) {			
			rcand.emplace_back(u);
		}
	}

	for (auto u : rcand) {
		removeFrCand(u);
	}
	//update excl
	for (ui i = 0; i < Excl.getSize(); i++) {
		ui u = Excl.get(i);
		if (!canMoveToP(u))
			rexcl.emplace_back(u);
	}
	for (auto u: rexcl) {
		Excl.remove(u);
	}
	branch(); // branch search	
	//recover
	PToCand(start);

	for (auto u : rcand) {
		addToCand(u);
	}
	for (auto u : rexcl) {
		Excl.add(u);
	}
}


void EnuBundle::checkSolution() {
	ui ismax;
	if (decompose)
		ismax = isGlobalMaximal();
	else
		ismax = 1;
	if (ismax) {
		cntplex++;
#ifdef SHOWSOL
		printf("Sol:");
		for (ui i = 0; i < P.getSize(); i++) {
			printf("%u ", bID[P.get(i)]);
		}
		printf("\n");
	/*
		ui rt = dbgCheckSolution();
		if (!rt)
			printf("Wrong solution pause\n");
			*/
#endif
	}
}

ui EnuBundle::dbgCheckSolution() {
	//check k-plex
	if (P.getSize() < lb) {
		printf("ERROR:Samller than lower bound\n");
		exit(-1);
	}
	//get sat
	set<ui> sat;
	for (ui i = 0; i < P.getSize(); i++) {
		ui u = P.get(i);
		ui deg = 0;
		for (ui j = 0; j < P.getSize(); j++) {
			if (P.get(j) != u && badc[u].test(P.get(j))) //trust block adjacency matrix
				deg++;
		}
		if (deg == P.getSize() - k)
			sat.insert(u);
		else if (deg < P.getSize() - k) {
			printf("ERROR:Not a k-plex\n");
			exit(0); //
		}
	}
	//check maximality
	for (ui u = 0; u < n; u++) {
		if (isInBlock(u) && P.contains(nID[u])) continue;// in subgraph
		if (core[u] + k >= P.getSize() + 1) {
			int scon = 0;
			int pcon = 0;
			for (ui j = pstart[u]; j < pstart[u + 1]; j++) {
				ui nei = edges[j];
				if (isInBlock(nei)) {
					ui idInBlk = nID[nei];
					if (P.contains(idInBlk))
						pcon++;
					if (sat.find(idInBlk) != sat.end())
						scon++;
				}
			}
			if (scon == sat.size() && pcon + k >= P.getSize() + 1) {
				printf("ERROR: Not maximal\n");
			}
				
		}
	}
	//check uniqueness
	vector<ui> vecsol;
	for (ui i = 0; i < P.getSize(); i++) {
		vecsol.emplace_back(bID[P.get(i)]);
	}
	Solution sol(vecsol);
	if (allsols.find(sol) != allsols.end()) {
		printf("ERROR: Not unique\n");
	}
	else {
		allsols.insert(sol);
	}
}

/**
Stop as the whole graph is a solution.
*/
void EnuBundle::stopAsSolution() { 
	vector<ui> rcand;
	vector<ui> rexcl;
	//nnodes++;
	for (ui i = 0; i < Cand.getSize(); i++)
		rcand.emplace_back(Cand.get(i));

	for (auto u : rcand) {
		CandToP(u);
	}
	for (ui i = 0; i < Excl.getSize(); i++) {
		ui u = Excl.get(i);
		if (!canMoveToP(u)) {
			rexcl.emplace_back(u);
		}
	}
#ifdef SHOWRECUR
	//printf("-------Nodes: %u (LEAF)-------------\n", nnodes);
	printf("P:");
	for (ui i = 0; i < P.getSize(); i++)
		printf("%u ", P.get(i));
	/*printf("\nCand:");
	for (ui i = 0; i < Cand.getSize(); i++)
		printf("%u ", Cand.get(i));
	printf("\nExcl:");
	for (ui i = 0; i < Excl.getSize(); i++)
		printf("%u ", Excl.get(i));
	printf("\n");*/
#endif
	//Check maximal, Excl is not changed
	if (rexcl.size() == Excl.getSize()) { // Excel will be empty
		checkSolution();
	}
	for (auto u : rcand) {
		PToCand(u);
	}
}

/*
ubr: suggestion a vertex for branch
*/
 void EnuBundle::branch() {
	 if (interrupt || (double)(clock() - sortclk) / CLOCKS_PER_SEC > maxsec) {
		 interrupt = 1;
		 return;
	 }
	nnodes++;
#ifdef SHOWRECUR
	printf("-------Nodes: %u-------------\n", nnodes);
	printf("P:");
	for (ui i = 0; i < P.getSize(); i++)
		printf("%u ", P.get(i));
	printf("\nCand:");
	for (ui i = 0; i < Cand.getSize(); i++)
		printf("%u ", Cand.get(i));
	printf("\nExcl:");
	for (ui i = 0; i < Excl.getSize(); i++)
		printf("%u ", Excl.get(i));
	printf("\n");
#endif
	if (Cand.empty() && Excl.empty()) {
		//maximal;
		//P.printList();
		if (P.getSize() >= lb) {
			checkSolution();
		}
		return;
	}
	if (Cand.empty()) {
#ifdef SHOWRECUR
		printf("|C|=0. Backtrack.\n");
		printf("--------BK-------------\n");
#endif
		return;
	}

	if (Cand.getSize() + P.getSize() < lb) {
#ifdef SHOWRECUR
		printf("|C|+|P|<lb. Backtrack.\n");
		printf("--------BK-------------\n");
#endif
		return;
	}

	//find mindeg 
	ui minu = UINT_MAX;
	/*
	for (ui u = 0; u < bn; u++) {			
		if (Cand.contains(u) || P.contains(u)) {
			if(minu == bn || neiInG[u] < neiInG[minu])
				minu = u;
		}
	}
	*/
	for (ui i = 0; i < P.getSize(); i++) {
		ui u = P.get(i);
		if (minu == UINT_MAX || neiInG[u] < neiInG[minu]) {
			minu = u;
		}
	}
	
	if (minu != UINT_MAX && neiInG[minu] + k < Cand.getSize() + P.getSize()) { // The whole graph can not be a k-plex
		ui szmax = k + neiInP[minu] - P.getSize();
		vector<ui> doing;
		for (ui i = 0; i < Cand.getSize(); i++) {
			ui u = Cand.get(i);
			if (u != minu && !badc[minu].test(u))
				doing.emplace_back(u);
		}		
		// assert(szmax < doing.size());
		multiRecurSearch(doing, szmax);
	}
	else {
		for (ui i = 0; i < Cand.getSize(); i++) {
			ui u = Cand.get(i);
			if (minu == UINT_MAX || neiInG[u] < neiInG[minu]) {
				minu = u;
			}
		}
		//Now, minu is the minimum degree vertex of G[P\cup C]
		if (neiInG[minu] >= Cand.getSize() + P.getSize() - k) {
			//The whole graph is a k-plex
			stopAsSolution();
#ifdef SHOWRECUR
			printf("G is a k-plex. Backtrack.\n");
			printf("--------BK-------------\n");
#endif
			return;
		}
		else if (Cand.getSize() + P.getSize() <= lb) {
#ifdef SHOWRECUR
			// the whole graph, which is of size lb, is not a k-plex,
			//thus it is not possible to find a better larger than lb
			printf("No better solution as G=lb, G is not kplex.Backtrack\n");
			printf("--------BK-------------\n");
#endif
			return;
		}		
		if (P.contains(minu)) {
			ui szmax = k + neiInP[minu] - P.getSize();
			vector<ui> doing;
			for (ui i = 0; i < Cand.getSize(); i++) {
				ui u = Cand.get(i);
				if (u != minu && !badc[minu].test(u))
					doing.emplace_back(u);
			}
			// assert(szmax < doing.size());
			multiRecurSearch(doing, szmax);
		}
		else { //minu in Cand
			//The first branch add minu to P
#ifdef SHOWRECUR
			printf("Add %u to P\n", minu);
#endif
			recurSearch(minu);

			//The second remove miu from Cand
			removeFrCand(minu);
			Excl.add(minu);
#ifdef SHOWRECUR
			printf("Remove %u\n", minu);
#endif
			branch(); 
			Excl.remove(minu);
			addToCand(minu);

		}
	}	
}
 
void EnuBundle::enumPlex(ui _k, ui _lb, uli _maxsec, ui _isdecompose, ui _quiete)
{	
	startclk = clock();
	k = _k;
	lb = _lb;
	maxsec = _maxsec;
	decompose = _isdecompose;
	quiete = _quiete;

	cntplex = 0;
	interrupt = 0; //interrupt the program when time out.
	
	Excl.init(n);
	P.init(n);
	Cand.init(n);
	neiInG = new ui[n];
	neiInP = new ui[n];
	
	dseq = new ui[n];
	core = new ui[n];
	dpos = new ui[n];	//vtxmp[v] is the position of v in degeneracy order 
	
	dist = new ui[n];
	common = new ui[n];
	szc1 = 0;
	cache1 = new ui[n];
	szc2 = 0;
	cache2 = new ui[n];
	szc3 = 0;
	cache3 = new ui[n];
	nID = new ui[n];
	bID = new ui[n];

	if (decompose) { // decompose the graph and enumerate k-plexes in smaller subgraphs
		ui maxcore = degeneracyOrder(dseq, core, dpos); // fast degeneracy order
		sortclk = clock();
		//build subgraph		
		ui maxbn = min(maxcore*maxcore + 1, n);
		bstart = new ui[maxbn + 1];
		bedges = new ui[m];
		badc = new MBitSet[maxbn];
		for (ui i = 0; i < maxbn; i++) {
			badc[i].allocacte(maxbn);
		}

		for (ui i = 0; i < n - 1; i++) {
			ui v = dseq[i];
			if (interrupt) break;
			if (core[v] + k >= lb) {
				ui sz = markBlock1(v, cache3); // don't use cache1 and cache3 as return adress
				if (sz < lb) {
					//printf("Vertex %u [%u] discard \n", i, v);
					continue;
				}
				else {
					buildBlock(v, cache3, sz);
				}
				if (!quiete)	printf("Block %u[%u]: %u %u\n", i, v, bn, bm);
#ifdef SHOWBLK	
				for (ui i = 0; i < bn; i++) {

					printf("%u(%u):", i, bID[i]);
					for (ui j = bstart[i]; j < bstart[i + 1]; j++) {
						printf("%u(%u) ", bedges[j], bID[bedges[j]]);
					}
					printf("\n");
				}
#endif
				//build cand;
				P.clear();
				//P.add(0);
				Cand.clear();

				memset(neiInP, 0, sizeof(int) * bn);
				for (ui u = 0; u < bn; u++) {
					Cand.add(u);
				}

				for (ui j = 0; j < bn; j++) {
					neiInG[j] = bstart[j + 1] - bstart[j];
				}
				Excl.clear();
				nnodes = 0;
				recurSearch(0);
				//branch();
				//Clear block
				if(!quiete) printf("Node number: %llu\n\n", nnodes);
			}
		}
	}
	else {
		bstart = new ui[n + 1];
		bedges = new ui[m];
		
		badc = new MBitSet[n];
		for (ui i = 0; i < n; i++) {
			badc[i].allocacte(n+1);
		}

		for (ui u = 0; u < n; u++) 
			cache3[szc3++] = u;
		buildBlock(0, cache3,szc3);

		sortclk = clock();
		printf("n=%u, n=%u\n", bn, bm);
		nnodes = 0;
		P.clear();
		Cand.clear();
		memset(neiInP, 0, sizeof(int) * bn);
		for (ui u = 0; u < bn; u++) {
			Cand.add(u);
		}

		for (ui j = 0; j < bn; j++) {
			neiInG[j] = bstart[j + 1] - bstart[j];
		}
		Excl.clear();
		branch();
	}
	enumclk = clock();
	printf("Number of %u-cplex larger than %u:  %u\n", k, lb, cntplex);
	printf("Total search time %.2f\n", Utility::elapse_seconds(startclk, enumclk));
	printf("Sort time %.2f\n", Utility::elapse_seconds(startclk, sortclk));
	//printf("Totoal nodes %u \n", nnodes);
}

EnuBundle::EnuBundle()
{
	n = 0;
	m = 0; 
	pstart = nullptr;
	edges = nullptr;
	reverse = nullptr;
	bstart = nullptr;
	bID = nullptr;
	bedges = nullptr;	
	badc = nullptr;
	binv = nullptr;	
}


EnuBundle::~EnuBundle()
{
	if (pstart != nullptr)
		delete[] pstart;
	if (edges != nullptr)
		delete[] edges;
	if (reverse != nullptr)
		delete[] reverse;
	delete[] dseq;
	delete[] dpos;
	delete[] core;
	delete[] nID;
	delete[] bstart;
	delete[] bedges;
	//delete[] badc;
	//delete binv;
	delete[] neiInG;
	delete[] neiInP;
	delete[] cache1;
	delete[] cache2;
	delete[] cache3;
}

