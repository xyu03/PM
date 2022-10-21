
#pragma once

#include <cstdio>
#include <cstdlib> 

/// <summary>
/// StaticBuckets: maps elements (0...maxid) to buckets (0...nbuckets)
/// Once in a bucket, an element cannot be moved or deleted.
/// And it's impossible to know which bucket an element belongs to.
/// </summary>
class StaticBuckets {
private:
    int maxid; //universe is 0...maxid
    int maxbucket; //buckets 0...maxbucket

    int *next;  //next[v]: successor of v in its list (-2:none, -1:end of list) 
    int *first; //first[b]: first element in b's bucket (-1: none)


public:
    void Reset() {
        int i;
        for (i = 0; i <= maxid; i++) next[i] = -2; //element is no list
        for (i = 0; i <= maxbucket; i++) first[i] = -1; //all buckets empty
    }

    void CheckEmpty() {
        fprintf (stderr, "Checking if StaticBuckets is empty... ");
        fprintf (stderr, "done.\n");
    }

    inline bool IsEmpty(int b) { return (first[b] == -1); }

    StaticBuckets (int _maxid, int _maxbucket) {
        maxid = _maxid;
        maxbucket = _maxbucket;

        next = new int[maxid + 1];
        first = new int[maxbucket + 1];

        Reset();
    }

	~StaticBuckets() {
		delete [] next;
		delete [] first;
	}

    void Reset(int b) {
        int i = first[b];
        if (i > -1) {
            do {
                int j = next[i];
                next[i] = -2;
                i = j;
            } while (i > 0);
        }
        first[b] = -1;
    }

    /// <summary>
    /// Insert an element into a bucket. Element must not be in any bucket already.
    /// </summary>
    /// <param name="e">the element</param>
    /// <param name="b">the bucket</param>
    void Insert(int e, int b) {
        next[e] = first[b];
        first[b] = e;
    }

    /// <summary>
    /// Enumerate all elements in bucket b.
    /// </summary>
    /// <param name="b">the bucket</param>
    /// <returns></returns>
    


	inline int GetFirst(int b) {return first[b];}
	inline int GetNext(int e) {return next[e];}
	
	
};

