package ru.ifmo.nds;

import ru.ifmo.nds.jfb.FenwickSweep;
import ru.ifmo.nds.jfb.RedBlackTreeSweep;
import ru.ifmo.nds.jfb.RedBlackTreeSweepHybridENS;
import ru.ifmo.nds.jfb.RedBlackTreeSweepHybridLinearNDS;

public class JensenFortinBuzdalov {
    private JensenFortinBuzdalov() {}

    private static final NonDominatedSortingFactory FENWICK_SWEEP = FenwickSweep::new;
    private static final NonDominatedSortingFactory RBTREE_SWEEP = RedBlackTreeSweep::new;
    private static final NonDominatedSortingFactory RBTREE_SWEEP_FNDS = RedBlackTreeSweepHybridLinearNDS::new;
    private static final NonDominatedSortingFactory RBTREE_SWEEP_ENS = RedBlackTreeSweepHybridENS::new;

    public static NonDominatedSortingFactory getRedBlackTreeSweepImplementation() {
        return RBTREE_SWEEP;
    }

    public static NonDominatedSortingFactory getFenwickSweepImplementation() {
        return FENWICK_SWEEP;
    }

    public static NonDominatedSortingFactory getRedBlackTreeSweepHybridFNDSImplementation() {
        return RBTREE_SWEEP_FNDS;
    }

    public static NonDominatedSortingFactory getRedBlackTreeSweepHybridENSImplementation() {
        return RBTREE_SWEEP_ENS;
    }
}
