package ru.ifmo.nds.jfb.hybrid;

import ru.ifmo.nds.jfb.HybridAlgorithmWrapper;
import ru.ifmo.nds.jfb.JFBBase;
import ru.ifmo.nds.util.ArrayHelper;
import ru.ifmo.nds.util.ArraySorter;
import ru.ifmo.nds.util.DominanceHelper;

import java.util.Arrays;

public class ENSH extends ENS {
    public ENSH(int threshold3D, int thresholdAll) {
        super(threshold3D, thresholdAll);
    }

    @Override
    public String getName() {
        return "ENSH (threshold 3D = " + threshold3D + ", threshold all = " + thresholdAll + ")";
    }

    @Override
    public HybridAlgorithmWrapper.Instance create(int[] ranks, int[] indices, double[][] points, double[][] transposedPoints) {
        return new Instance(ranks, indices, points, threshold3D, thresholdAll);
    }

    private static final class Instance extends ENS.Instance {
        protected Instance(int[] ranks, int[] indices, double[][] points, int threshold3D, int thresholdAll) {
            super(ranks, indices, points, threshold3D, thresholdAll);
        }

        private int transplantRanksAndCheckWhetherAllAreSame(int goodFrom, int goodUntil, int ranksAndSlicesOffset, int sortedIndicesOffset) {
            int curI = sortedIndicesOffset;

            int min = ranks[indices[goodFrom]];
            for (int i = goodFrom; i < goodUntil; ++i) {
                int rank = ranks[indices[i]];
                if (min > rank) {
                    min = rank;
                    curI = sortedIndicesOffset;
                }
                if (min == rank) {
                    space[curI++] = indices[i];
                }
            }
            space[ranksAndSlicesOffset++] = min;

            for (int i = goodFrom, ri = ranksAndSlicesOffset, si = curI; i < goodUntil; ++i) {
                int rank = ranks[indices[i]];
                if (rank == min)
                    continue;
                space[si++] = ri;
                space[ri++] = -rank;
            }
            return curI;
        }

        @Override
        public int helperBHook(int goodFrom, int goodUntil, int weakFrom, int weakUntil, int obj, int tempFrom, int maximalMeaningfulRank) {
            int goodSize = goodUntil - goodFrom;
            if (notHookCondition(goodSize + weakUntil - weakFrom, obj)) {
                return -1;
            }

            int sortedIndicesFrom = tempFrom * STORAGE_MULTIPLE;
            int sortedIndicesUntil = sortedIndicesFrom + goodSize;

            int ranksAndSlicesOffset = sortedIndicesUntil;
            int sliceOffset = ranksAndSlicesOffset + goodSize;

            int minRankFrom = sortedIndicesFrom;
            int minRankCurUntil = sortedIndicesFrom;
            int minRankUntil = sortedIndicesFrom = transplantRanksAndCheckWhetherAllAreSame(goodFrom, goodUntil, ranksAndSlicesOffset, sortedIndicesFrom);
            int minRank = space[ranksAndSlicesOffset++];

            if (sortedIndicesFrom == sortedIndicesUntil) {
                return helperBSingleRank(minRank, goodFrom, goodUntil, weakFrom, weakUntil, obj, maximalMeaningfulRank, tempFrom);
            }

            ArraySorter.sortIndicesByValues(space, space, sortedIndicesFrom, sortedIndicesUntil);
            int sliceLast = distributePointsBetweenSlices(space, sortedIndicesFrom, sortedIndicesUntil, sliceOffset, tempFrom);
            int minOverflowed = weakUntil;
            for (int weak = weakFrom, good = goodFrom, sliceOfGood = ranksAndSlicesOffset; weak < weakUntil; ++weak) {
                int wi = indices[weak];
                for (; good < goodUntil && indices[good] < wi; ++good) {
                    if (ranks[indices[good]] == minRank)
                        continue;
                    exPoints[space[space[sliceOfGood++] + 1]++] = points[indices[good]];
                }
                int weakRank = findRankInSlices(sliceOffset, sliceLast, wi, obj, sortedIndicesFrom);
                if (minRank >= weakRank) {
                    minRankCurUntil = ArrayHelper.findWhereNotSmaller(space, minRankCurUntil, minRankUntil, wi);
                    double[] point = points[wi];
                    for (int i = minRankFrom; i < minRankCurUntil; ++i) {
                        if (DominanceHelper.strictlyDominatesAssumingLexicographicallySmaller(points[space[i]], point, obj)) {
                            ranks[wi] = weakRank = minRank + 1;
                            break;
                        }
                    }
                }
                if (weakRank > maximalMeaningfulRank && minOverflowed > weak) {
                    minOverflowed = weak;
                }
            }
            Arrays.fill(exPoints, tempFrom, tempFrom + goodUntil - goodFrom, null);
            return JFBBase.kickOutOverflowedRanks(indices, ranks, maximalMeaningfulRank, minOverflowed, weakUntil);
        }
    }
}
