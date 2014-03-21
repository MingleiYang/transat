/*
 * Copyright: Nicholas P Wiebe (2008-2009) and Irmtraud M Meyer (2008-2009)
 * License: licensed under the GNU General Public License version 3 (GPLv3)
 * HelixFinder.h
 *
 * Finds and stores helices in a given alignment
 */

#ifndef HELIXFINDER_H_
#define HELIXFINDER_H_

#include <vector>
#include <list>
#include "AlignedHelix.h"
#include "Tree.h"

using namespace std;

class Alignment;

enum Realigner { NO_REALIGN, TCOFFEE };

class HelixFinder {
public:
	HelixFinder(Alignment * alignment_);
	virtual ~HelixFinder();

	void findAllHelices();

	void findAllHelicesGrow(Tree & tree, double seedThreshold, double growThreshold);

	void findCompetingHelices();

	/**
	 * not implemented
	 */
	void findTrueHelices();

	/**
	 * idea: if two helices overlap, we would like to throw away the lowest of these helices
	 * actions on a sorted list:
	 * highest dominates all others which overlap at all
	 * would like to shrink others rather than throw them away
	 */

	/**
	 * calculate pvalues for log ratios of all found helices,
	 * by generating randomized alignments to simulate null distribution.
	 *
	 * @param randomSamples Number of randomized alignments to generate
	 * @param doPvalues if false, don't do pvalue calculation - just fill that column with zeros
	 */
	void allHelicesPvalueTable(int randomSamples, Tree & tree, bool doPvalues=true);

	void balancedSparseHelixTable(Tree & tree);

	void sparseHelixTable(Tree & tree);

	/**
	 * calculates the coverage of the consensus structure (i.e. the fraction of bps in the consensus structure that
	 * can be mappped to at least one found helix)
	 *
	 * @param trueHelixCutoff A helix must be composed of greater than this fraction of consensus bases to be considered a match for those bases.
	 * if trueHelixCutoff is set to 1.0, then only helices that are entirely composed of consensus bases are considered
	 * @param coverages Sets this to the fraction of bps in the consensus structure that are covered by helices that meet the trueHelix cutoff criterion
	 * @param exactMatchCoverage Sets this to the fraction of bps in the consensus structre that are covered by exact matches (See AlignedHelix::isConsensusHelix)
	 */
	void coverage(double & coverage, double & exactMatchCoverage, double trueHelixCutoff = 0.5);

	vector<AlignedHelix> helices;

	static Realigner realign;
	static bool verbose_out;


private:
	Alignment * alignment; //reference instead of pointer?

	/**
	 * list of
	 */
	list<pair<AlignedHelix*, AlignedHelix*> > competesWith;

	/**
	 * find all helices in a specific sequence from the alignment
	 * (not added to helices vector, just to helicesByOuterBP vector)
	 */
	void findAllHelices(int seqIndex, vector<list<AlignedHelix> > & helicesByOuterBP);

	/**
	 * converts a pair of positions to an index in a flattened UT matrix, as used in helicesByOuterBP
	 */
	int convert(int pos5, int pos3);

	void clearDuplicates(vector<list<AlignedHelix> > & helicesByOuterBP);
};

#endif /* HELIXFINDER_H_ */
