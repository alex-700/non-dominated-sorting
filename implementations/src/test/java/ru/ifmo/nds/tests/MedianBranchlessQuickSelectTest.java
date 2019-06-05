package ru.ifmo.nds.tests;

import ru.ifmo.nds.util.ArrayHelper;

public class MedianBranchlessQuickSelectTest extends MedianTestsBase {
    @Override
    protected double destructiveMedian(double[] array, int until) {
        return ArrayHelper.destructiveBranchlessMedian(array, 0, until, new double[array.length]);
    }
}
