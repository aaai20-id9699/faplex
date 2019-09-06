#pragma once
#include "Defines.h"
#include <algorithm>
#include <cassert>
class RandList
{
private:
	ui *vlist;
	ui *vpos;
	ui vnum;
	ui cap;
public:
	RandList() {
		vlist = vpos = nullptr;
		//vnum = cap = 0;
	};
	void init(int _cap) {		
		cap = _cap;
		if (vlist == nullptr) 
			vlist = new ui[cap];
		if (vpos == nullptr) 
			vpos = new ui[cap];

		vpos = new ui[cap];
		vnum = 0;
		for (ui i = 0; i < cap; i++) {
			vpos[i] = cap;
		}
	}	
	void add(int vid) {
		assert(vpos[vid] == cap);
		vlist[vnum] = vid;
		vpos[vid] = vnum;
		vnum++;

	};
	void remove(int vid) {
		assert(vpos[vid] < vnum);
		ui last_id = vlist[vnum - 1];
		ui id_pos = vpos[vid];
		vlist[id_pos] = last_id;
		vpos[last_id] = id_pos;
		vnum--;
		vpos[vid] = cap; /*set as visited*/
	}
	void clear() {			
		for(ui i = 0; i < vnum; i++)
			vpos[i] = cap;
		vnum = 0;
	}
	ui get(ui i) {
		assert(i < vnum);
		return vlist[i];
	}
	bool contains(int vid) {
		return vpos[vid] != cap;
	}
	bool empty() { return vnum == 0; }
	ui getSize() { return vnum; }
	ui getCap() { return cap; }
	void dispose() {
		if (vlist != nullptr) {
			delete[] vlist;
			vlist = nullptr;
		}
		if (vpos != nullptr) {
			delete[] vpos;
			vpos = nullptr;
		}
	}
	~RandList() {
		dispose();
	}
#ifdef DBGMOD		
	void printList(FILE *f = stdout) {
		fprintf(f, "Total %d: ", vnum);
		int *tmp_lst = new int[cap];
		memcpy(tmp_lst, vlist, vnum * sizeof(int));
		std::sort(tmp_lst, tmp_lst + vnum);
		//qsort(tmp_lst, vnum, sizeof(int), cmpfunc);
		for (ui i = 0; i < vnum; i++) {
			fprintf(f, "%d ", tmp_lst[i]);
		}
		fprintf(f, "\n");
	};
#else
	void printList(FILE *f = stdout) {};
#endif
};

extern int *twoPow;

class MBitSet{
private:
	int n, m;
	int cap;
	unsigned *buf;
	
	const int msk = (1 << 16) - 1;	
	
	// int twoPow[1 << 16];
	int lb(unsigned x) {
		if (x & msk)
			return twoPow[x & msk];
		return twoPow[x >> 16 & msk] + 16;
	}
	
public:
	MBitSet() {		
		buf = nullptr;
		cap =  n = m = 0;
	}
	~MBitSet() {
		if (buf != nullptr)
			delete[] buf;
	}
	void allocacte(int _cap) {
		cap = _cap;

		buf = new unsigned[(cap >>5) + 1];
	}
	void reinit(int _n) {
		assert(_n <= cap);
		m = _n & 31;
		n = _n >> 5;		
		//buf = new unsigned[n + 1];
		for (int i = 0; i <= n; ++i)
			buf[i] = 0;		
		//for (int i = 0; i < 16; ++i)
		//	twoPow[1 << i] = i;//Ԥ�����2^i 
	}
	//FLIP all the bits
	void flip() {
		for (int i = 0; i < n; ++i)
			buf[i] = ~buf[i];
		buf[n] ^= ((unsigned)1 << m) - 1;
	}
	void set(int x) {
		//assert(x < cap);
		buf[x >> 5] ^= (unsigned)1 << (x & 31);
	}
	bool test(int x) {
		return buf[x >> 5] >> (x & 31) & 1;
	}
	
	//return the lowest non-zero bit
	int lowbit() {
		for (int i = 0; i <= n; ++i) {
			unsigned x = buf[i] & ((~buf[i]) + 1);
			if (x)
				return lb(x) + (i << 5);
		}
		return -1;
	}
	bool empty() {
		for (int i = 0; i <= n; ++i)
			if (buf[i])
				return false;
		return true;
	}
	void operator &=(const MBitSet &rhs) {
		for (int i = 0; i <= n; ++i)
			this->buf[i] &= rhs.buf[i];
	}
	//copy the bitset INTO rhs
	void copy(MBitSet &rhs) const {
		rhs.n = n;
		rhs.m = m;
		if (rhs.buf == NULL)
			rhs.buf = new unsigned[n + 1];
		for (int i = 0; i <= n; ++i)
			rhs.buf[i] = buf[i];
	}

};

