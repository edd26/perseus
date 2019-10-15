/*
 * Morse.h
 * Routines to perform Discrete Morse Theoretic reductions via free coface
 * collapses of a Cell Complex
 */

#ifndef MORSE_H_
#define MORSE_H_

#include <string>
#include <deque>


#define debout cout

#include "Complex.h"

// morse complex
template <typename C, typename BT = int>
class MComplex: public Complex<C,BT>
{
public:

	// critical cell info
	// the index of the next potential candidate cell for being critical
	num minindex;
	bool isred; // are we reducing or coreducing?
	bool hyperq; // are we being aggressive with the queues?

	MComplex():Complex<C,BT>()
	{
	    hyperq = false;
	    isred = false;
		minindex = 0;
	}

	~MComplex()
	{
		//destroyPersData();
	}

	// cell insertion routines
	void quicksert(const vector<Cell<C,BT>*>&); // inserts a vector of cells into the complex
	bool insertCell(Cell<C,BT>*); // inserts a single cell into complex

	pair<num,bool> getMinDim(const BT&); // re-evaluates minimum uncritical dimension for a given birth time
    pair<num,bool> getMaxDim(const BT&); // re-evaluates maximum uncritical dimension for a given birth time

	// cell removal routines
	bool removeCell(Cell<C,BT>* tokill, typename CLIST::iterator& pos, bool asCrit=false); // removes a single cell from complex

	// generator output!
	void showGens (const BT& born, const num& dim, bool resbirth = false, ostream& out = cout) const; // output generators

	void moveOver_Sort (MComplex<C,BT>&,bool);
	void moveOver (MComplex<C,BT>&);
	void moveOver_Shuffle (MComplex<C,BT>&);

	double boundaryMatrix (const num& updim, vector<vector<C> >& bmat, BT birth = BANBT)const;
    double getBoundaryDensity (const num& updim, BT birth = BANBT) const;
    double getUnitFraction (const num& updim, BT birth = BANBT) const;


	// **************** COREDUCTION ROUTINES, SEE ALGOS/MORSEREDUCTION.HPP ***********************
	pair<num,int> CoReduce(const BT& birth = INITBT, const BT& ignore = BANBT,
                           bool deep = true, bool tracer = false);	// coreduces given frame iteratively with
                                                        // lowest dim cell as seed
	pair<num,bool> CoredMorsify(const BT& birth, bool deep = true,
                                BT except = BANBT,bool tracer=false); // reduces complex via iterated coreduction, removes uncritical cells from
	                                                // memory if deep is on.
	num MorseCoreduce(bool savegens = false, bool tracer = false); // coredmorsify + ff collapses, top level morse reduction
	num MorseCoreduceExcept(const BT& start = INITBT, const BT& exception = BANBT, bool savegens = false); // same as above, but leaves given birth frame unreduced


	void makeCritical_Cored(Cell<C,BT>* tocrit, bool makekgen = false); // processes given cell as critical i.e. makes an Ace...
	bool markPair_Cored(Cell<C,BT>* king, Cell<C,BT>* queen, deque<Cell<C,BT>*>& Queue, num rootdim,
                        typename CLIST::iterator& kpos, typename CLIST::iterator& qpos, bool makekgens=false); // processes KQ pair
	bool markUncritical_Cored(const BT&); // marks all cells uncritical for repeating coredmorsify
	// free face collapse routines
	pair<num,int> freeFaceCollapse(const BT&,Cell<C,BT>*&,bool);
	void removeFreeFace(Cell<C,BT>* torem, bool deep=true);
	num getNextCrit_Cored(const BT& birth, Cell<C,BT>*& newace, bool deep = true);

	// gradient path and generator handling functions
	void updateGpathNearAce_Cored(Cell<C,BT>*);
	void updateGpathNearQueen_Cored(Cell<C,BT>*, deque<Cell<C,BT>*>&, bool);
	bool inheritGpathFromKing_Cored(Cell<C,BT>*, Cell<C,BT>*,bool);
	void updateKgenNearQueen_Cored(Cell<C,BT>*, Cell<C,BT>*, const C&);



    // **************** REDUCTION ROUTINES, SEE ALGOS/MORSEREDUCTION.HPP ***********************

    pair<num,int> Reduce(const BT& birth = INITBT, const BT& ignore = BANBT,
                         bool deep = true, bool tracer = false);	// reduces given frame iteratively with highest
                                                        // dim cell as seed
    pair<num,bool> RedMorsify(const BT& birth, BT except, bool tracer, bool deep = true); // reduces complex via iterated coreduction, removes uncritical cells from
	                                                // memory if deep is on.
	num MorseReduce(bool savegens = false, bool tracer = false); // redmorsify + ff collapses, top level morse reduction
	num MorseReduceExcept(const BT&,const BT&,bool); // same as above, but leaves given birth frame unreduced


    void makeCritical_Red(Cell<C,BT>* tocrit, bool makekgen = false); // processes given cell as critical i.e. makes an Ace...
	bool markPair_Red(Cell<C,BT>* queen, Cell<C,BT>* king, deque<Cell<C,BT>*>& Queue, num rootdim,
                      typename CLIST::iterator& qpos, typename CLIST::iterator& kpos, bool makekgens=false); // processes QK pair
	bool markUncritical_Red(const BT&); // marks all cells uncritical for repeating redmorsify
    void removeFreeCoface(Cell<C,BT>* torem, bool deep=true);

    num getNextCrit_Red(const BT& birth, Cell<C,BT>*& newace, bool deep = true); // next critical cell for reduction

    // gradient path and generator handling functions
	void updateGpathNearAce_Red(Cell<C,BT>*);
	void updateGpathNearKing_Red(Cell<C,BT>*, deque<Cell<C,BT>*>&, bool);
	bool inheritGpathFromQueen_Red(Cell<C,BT>*, Cell<C,BT>*, bool);
	void updateKgenNearKing_Red(Cell<C,BT>*, Cell<C,BT>*, const C&);



    // ********************common to reduction and coreduction**********************
    bool setMinIndex(const BT);
    void enqueue (const CCHAIN&, deque<Cell<C,BT>*>&);
	void memorySaver(BT birth = INITBT);
	bool verify(const BT&);
	bool hasUncritCells(const BT&, const num&) const;
	bool hasUncritCells(const BT&) const;
	num freeFaceCollapser();
	num freeCofaceCollapser();
	num Alternator(const BT& birth = INITBT, bool savegens = false);

    // ******************* high level wrappers****************************

    num TopAlternator(bool redfirst=true, double thresh = 0.2,
                       bool savegens = false, bool tracer = false);
    num HybReduce(bool&,bool,const BT); // hybrid reduce/coreduce strategy
    pair<num,bool> HybMorsify(bool&, const BT&, bool, const BT);
	void MorseWrapper_Hybrid(bool,double); // hybrid strategy
	num MorseWrapper_Cored(bool savegens = false, double thresh = 0.3); // iterated morse coreducer!
	void MorseWrapperExcept_Cored(const BT& start, const BT& except, bool savegens = false,
                                  double thresh = 0.3);
	num MorseWrapper_Red(bool savegens = false, double thresh = 0.3); // iterated morse coreducer!
	void MorseWrapperExcept_Red(const BT& start, const BT& except, bool savegens = false,
                                double thresh = 0.3);
};

#endif /* MORSE_H_ */
