
#pragma once

#include <cstdio>
#include <cstdlib>


class UniverseSet {
	private:
		int *perm;
		int *pos; //pos[v]: position of element v
		int maxid;  //maximum allowed identifier
		int nextpos; //position of next insertion
		int nil;

    
		/// <summary>
		/// Swap element in position pi into position nextpos.
		/// </summary>
		/// <param name="pi">current element position</param>
		void MakeNext(int pi) {
			int i = perm[pi]; //element originally in position i
			int j = perm[nextpos]; //element originally in nextpos

			//i goes to nextpos
			perm[nextpos] = i;
			pos[i] = nextpos;

			//j goes to pi
			perm[pi] = j;
			pos[j] = pi;
		}

		void die (const string &msg) {
			fprintf(stderr, "(UniverseSet): %s.\n", msg.c_str());
		}


	public:
		/// <summary>
		/// Create a set representing elements 0...n.
		/// </summary>
		/// <param name="n"></param>
		UniverseSet (int _maxid) {
			maxid = _maxid;
			nil = maxid + 1;
			perm = new int [_maxid+1];
			pos = new int [_maxid+1];
			HardReset();
		}

		inline int PickAny() {
			if (IsEmpty()) die ("Cannot pick element from an empty set.");
			return (perm[0]);
		}

		/// <summary>
		/// Empty the set.
		/// </summary>
		void HardReset()
		{
			for (int i=0; i<=maxid; i++) {
				pos[i] = i;
				perm[i] = i;
			}
			nextpos = 0;
		}

		void Reset() {
			while (!IsEmpty()) {Remove(perm[0]);}
		}


		inline bool Contains(int i) const {return (pos[i]<nextpos);}
		inline bool IsEmpty() {return (nextpos==0);}

		inline bool Insert(int i) {
			int pi = pos[i];
			if (pi < nextpos) return false; //already there: nothing to do
			MakeNext(pi);
			nextpos++;
			return true;
		}

		inline bool Remove (int i) {
			int pi = pos[i];
			if (pi >= nextpos) return false;
			nextpos--;
			MakeNext(pi);
			return true;
		}



    inline int Count() {return (nextpos);}

    inline int MaxId() {return maxid;}

    inline void Copy(UniverseSet *s) {
        if (s->MaxId() != MaxId()) die("Cannot copy from incompatible set");

        Reset(); //make this set empty
		Unite(s); //take the union with the other set
    }

	inline void Unite (UniverseSet *s) {
		if (s->MaxId() > MaxId()) die ("Cannot unite potentially larger set.\n");
		int i, end;
		for (s->GetBounds(i,end); i<end; i++) {
			Insert(s->PickPos(i));
		}

	}

	inline void GetBounds(int &start, int &end) {
		start = 0;
		end = nextpos;
	}

	inline int PickPos(int p) {
		if (p<0 || p>nextpos) die ("picking from invalid position");
		return perm[p];
	}



    void Output() {
        int i, end;
        for (i=0; i<=maxid; i++) {
            if (Contains(i)) fprintf (stderr, "#");
            else fprintf (stderr, "_");
        }
		for (GetBounds(i,end); i<end; i++) {
			fprintf (stderr, " %d", PickPos(i));
		}
		fprintf (stderr, "\n");

        
        CheckConsistency();
    }

    void CheckConsistency() {
        fprintf (stderr, "Check consistency.");

        for (int i=0; i<=maxid; i++) {
            if (perm[pos[i]] != i) fprintf(stderr, "Inconsistency at position %d.\n", i);
            if (pos[perm[i]] != i) fprintf(stderr, "Inconsistency for element %d.\n", i);
        }
    }

	~UniverseSet() {
		delete [] perm;
		delete [] pos;
	}
};

