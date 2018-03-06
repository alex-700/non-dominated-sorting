package ru.ifmo.nds.jfb;

import ru.ifmo.nds.bos.ImprovedAdaptedForHybrid;

import java.util.Arrays;
import java.util.Random; // TODO delete

public class RedBlackTreeSweepHybridBOS extends RedBlackTreeSweep {
    private final ImprovedAdaptedForHybrid bos;
    private double[][] tempPoints;
    private int[] tempRanks;
    private Random random = new Random(239); // TODO delete

    private static final int THRESHOLD_3D = 20;
    private static final int THRESHOLD_ALL = 200;

    public RedBlackTreeSweepHybridBOS(int maximumPoints, int maximumDimension, int allowedThreads) {
        super(maximumPoints, maximumDimension, allowedThreads);
        bos = new ImprovedAdaptedForHybrid(maximumPoints, maximumDimension);
        tempPoints = new double[maximumPoints][maximumDimension];
        tempRanks = new int[maximumPoints];
    }

    @Override
    public String getName() {
        return "Jensen-Fortin-Buzdalov sorting, "
                + getThreadDescription()
                + " (tree sweep, hybrid with Best Order Sort)";
    }

    @Override
    protected boolean helperAHookCondition(int size, int obj) {
//        switch (obj) {
//            case 1: return false;
//            case 2: return size < THRESHOLD_3D;
//            default: return size < THRESHOLD_ALL;
//        }
        return false;
    }

    @Override
    protected int helperAHook(int from, int until, int obj) {
        getPoints(from, until, obj + 1, tempPoints);
        getRanks(from, until, tempRanks);

        bos.sortCheckedWithRespectToRanks(
                tempPoints,
                tempRanks,
                until - from,
                obj + 1,
                maximalMeaningfulRank);

        for (int i = from; i < until; i++) {
            ranks[indices[i]] = tempRanks[i - from];
        }
        return until;
    }

    @Override
    protected boolean helperBHookCondition(int goodFrom, int goodUntil, int weakFrom, int weakUntil, int obj) {
//        return random.nextBoolean(); // TODO fix некоторые тесты парают !
        return goodFrom != 0;
//        return true;
    }

    @Override
    protected int helperBHook(int goodFrom, int goodUntil, int weakFrom, int weakUntil, int obj, int tempFrom) {
        for(int i = 0; i < getMaximumPoints(); i++) { // TODO delete
            Arrays.fill(tempPoints[i], -1);
        }
        Arrays.fill(tempRanks, -1); // TODO delete
        getPoints(0, Math.max(goodUntil, weakUntil), obj + 1, tempPoints); // TODO fix max
        getRanks(0, Math.max(goodUntil, weakUntil), tempRanks); // TODO fix max

        int newWeakUntil = bos.sortCheckedWithRespectToRanksHelperB(
                tempPoints,
                tempRanks,
                goodFrom,
                goodUntil,
                weakFrom,
                weakUntil,
                obj + 1,
                maximalMeaningfulRank);

        for (int i = goodFrom; i < goodUntil; i++) {
            ranks[indices[i]] = tempRanks[i];
        }

        for (int i = weakFrom; i < weakUntil; i++) {
            ranks[indices[i]] = tempRanks[i];
        }

        if(newWeakUntil == -1) {
            throw new IllegalStateException("Invalid newWeakUntil");
        }
        return newWeakUntil;
    }
}

