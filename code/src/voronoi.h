
#pragma once

/*----------------------------------------------------
 | Information necessary to maintain Voronoi diagrams
 *---------------------------------------------------*/

struct VoronoiUnit {
public:
	int vorbase;   //base of voronoi region
    int parc;      //label of the parent edge (WARNING: CHANGE NAME TO PEDGE)
    EdgeCost dist; //distance from the base

	inline void Reset() { 
		vorbase = parc = 0;
		dist = 0; 
	}

    inline void Reset(EdgeCost d) {
		vorbase = parc = 0;
        dist = d;
	}

    inline void Update(int b, int p, EdgeCost d) {
		vorbase = b;
        parc = p;
        dist = d;
	}

    inline void CopyFrom(VoronoiUnit x) {
		vorbase = x.vorbase;
        parc = x.parc;
        dist = x.dist;
    }
};


class VoronoiData {
private:
	VoronoiUnit *vor;
    int n;

public:
    void Reset() {
		for (int v = 0; v <= n; v++) {vor[v].Reset();}
	}

    void CopyFrom(int v, VoronoiData &other) {
		Update(v, other.GetBase(v), other.GetParentArc(v), other.GetDistance(v));
	}

    void Clone(VoronoiData &other) {
		for (int v=0; v<=n; v++) CopyFrom(v, other);
	}

    inline void Reset(int v) {vor[v].Reset();}

    VoronoiData(int _n) {
		n = _n;
        vor = new VoronoiUnit[n + 1];
        Reset();
    }

	~VoronoiData() {delete [] vor;}

    inline void Update(int v, int b, int p, EdgeCost d) {vor[v].Update(b, p, d);}

	// create voronoi region for v
    inline void MakeBase(int v) {vor[v].Update(v,0,0);}

	// query functions
	inline int GetBase(int v) const {return vor[v].vorbase;}
    inline EdgeCost GetDistance(int v) const {return vor[v].dist;}
	inline int GetParentArc(int v) const {return vor[v].parc;}
};

