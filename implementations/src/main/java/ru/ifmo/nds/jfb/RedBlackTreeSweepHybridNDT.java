package ru.ifmo.nds.jfb;

import ru.ifmo.nds.ndt.Split;
import ru.ifmo.nds.ndt.SplitBuilder;
import ru.ifmo.nds.ndt.TreeRankNode;

public class RedBlackTreeSweepHybridNDT extends RedBlackTreeSweep {
    private static final int THRESHOLD_3D = 100;
    private static final int THRESHOLD_ALL = 20000;

    private SplitBuilder splitBuilder;
    private TreeRankNode tree;
    private final int threshold;

    public RedBlackTreeSweepHybridNDT(int maximumPoints, int maximumDimension, int threshold) {
        super(maximumPoints, maximumDimension, 1);

        this.threshold = threshold;
        this.splitBuilder = new SplitBuilder(maximumPoints);
        this.tree = TreeRankNode.EMPTY;
    }

    @Override
    protected void closeImpl() {
        super.closeImpl();
        splitBuilder = null;
        tree = null;
    }

    @Override
    public String getName() {
        return "Jensen-Fortin-Buzdalov sorting, " + getThreadDescription() + " (tree sweep, hybrid with ENS-NDT, threshold " + threshold + ")";
    }

    @Override
    protected boolean helperAHookCondition(int size, int obj) {
        switch (obj) {
            case 1: return false;
            case 2: return size < THRESHOLD_3D;
            default: return size < THRESHOLD_ALL;
        }
    }

    @Override
    protected boolean helperBHookCondition(int goodFrom, int goodUntil, int weakFrom, int weakUntil, int obj) {
        return helperAHookCondition(goodUntil - goodFrom + weakUntil - weakFrom, obj);
    }

    @Override
    protected int helperAHook(int from, int until, int obj) {
        Split split = splitBuilder.result(transposedPoints, from, until, indices, obj + 1, threshold);

        int minOverflow = until;
        int M = obj + 1;
        tree = TreeRankNode.EMPTY;
        for (int i = from; i < until; ++i) {
            int idx = indices[i];
            ranks[idx] = tree.evaluateRank(points[idx], ranks[idx], split, M);

            if (ranks[idx] <= maximalMeaningfulRank) {
                tree = tree.add(points[idx], ranks[idx], split, threshold);
            } else if (minOverflow > i) {
                minOverflow = i;
            }
        }
        tree = null;
        return kickOutOverflowedRanks(minOverflow, until);
    }

    @Override
    protected int helperBHook(int goodFrom, int goodUntil, int weakFrom, int weakUntil, int obj, int tempFrom) {
        Split split = splitBuilder.result(transposedPoints, goodFrom, goodUntil, indices, obj + 1, threshold);

        int minOverflow = weakUntil;
        int M = obj + 1;
        tree = TreeRankNode.EMPTY;
        for (int good = goodFrom, weak = weakFrom; weak < weakUntil; ++weak) {
            int wi = indices[weak];
            int gi;
            while (good < goodUntil && (gi = indices[good]) < wi) {
                tree = tree.add(points[gi], ranks[gi], split, threshold);
                ++good;
            }
            ranks[wi] = tree.evaluateRank(points[wi], ranks[wi], split, M);
            if (minOverflow > weak && ranks[wi] > maximalMeaningfulRank) {
                minOverflow = weak;
            }
        }
        tree = null;
        return kickOutOverflowedRanks(minOverflow, weakUntil);
    }
}
