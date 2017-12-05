package scratch.kevin.griddedInterpGMPE;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIMR;

import com.google.common.base.Preconditions;

public abstract class AbstractGMPEInterpolation<E> implements Iterable<E> {
	
	private final String name;
	
	private static final double PRECISION_SCALE = 1 + 1e-14;
	
	public static abstract class EvenlySpacedDouble extends AbstractGMPEInterpolation<Double> {
		
		private final double min, max, delta;
		private final int numBins;
		private final boolean logGridding;
		private final boolean logXInterp;

		EvenlySpacedDouble(String name, double min, double max, int numBins,
				boolean logXInterp, boolean logGridding) {
			super(name);
			if (logGridding) {
				this.min = Math.log(min);
				this.max = Math.log(max);
			} else {
				this.min = min;
				this.max = max;
			}
			this.numBins = numBins;
			Preconditions.checkArgument(min <= max, "Min must be <= max");
			if (min < max) {
				this.delta = (max - min) / (numBins - 1);
			} else {
				// min == max
				Preconditions.checkArgument(numBins == 1, "Num bins must be 1 if min == max");
				this.delta = 0;
			}
			this.logGridding = logGridding;
			this.logXInterp = logXInterp;
		}

		@Override
		public int getNumBins() {
			return numBins;
		}

		@Override
		public Double getValue(int index) {
			if (index < 0 || index > ( numBins -1 ))
				throw new IndexOutOfBoundsException("no point at index "+index);
			if (logGridding)
				return Math.exp((min + delta*index));
			return min + delta*index;
		}
		
		private int getClosestIndex(double value) {
			if (logGridding)
				value = Math.log(value);
			double iVal = PRECISION_SCALE * (value - min) / delta;
			int i = (delta == 0) ? 0 : (int) Math.round(iVal);
			return (i<0) ? 0 : (i>=numBins) ? numBins-1 : i;
		}
		
		private int getIndexBefore(double value) {
			int index = getClosestIndex(value);
			if (getValue(index) > value)
				index--;
			return index;
		}

		@Override
		public double getInterpolatedBinIndex(Double value) {
			int binBefore = getIndexBefore(value);
			if (binBefore < 0)
				return 0;
			if (binBefore == numBins - 1)
				return numBins - 1;
			double valueBefore = getValue(binBefore);
			if (value == valueBefore)
				return binBefore;
			double valueAfter = getValue(binBefore+1);
			// actual interp
			double beforeDelta;
			double afterDelta;
			if (logXInterp) {
				beforeDelta = Math.log(value) - Math.log(valueBefore);
				afterDelta = Math.log(valueAfter) - Math.log(value);
			} else {
				beforeDelta = value - valueBefore;
				afterDelta = valueAfter - value;
			}
			Preconditions.checkState(beforeDelta >= 0, "Bad beforeDelta=%s for %s with range [%s %s], log=%s, binBefore=%s",
					beforeDelta, value, min, max, logXInterp, binBefore);
			Preconditions.checkState(afterDelta >= 0, "Bad afterDelta=%s for %s with range [%s %s], log=%s, binBefore=%s",
					afterDelta, value, min, max, logXInterp, binBefore);
//			System.out.println(value+": "+binBefore+" "+beforeDelta+" "+afterDelta);
			double index = binBefore + beforeDelta/(beforeDelta + afterDelta);
			Preconditions.checkState(index >= 0 && index <= numBins-1, "Index outside bounds. value=%s, binBefore=%s, beforeDelta=%s",
					value, binBefore, beforeDelta);
//			Preconditions.checkState(getValue((int)index) <= value);
//			if (index < numBins-1)
//				Preconditions.checkState(getValue((int)index+1) >= value);
			return index;
		}
		
		public double getMin() {
			return min;
		}
		
		public double getMax() {
			return max;
		}
	}
	
	public static abstract class Discrete<E> extends AbstractGMPEInterpolation<E> {
		
		private List<E> values;

		Discrete(String name, List<E> values) {
			super(name);
			
			this.values = Collections.unmodifiableList(values);
		}

		@Override
		public int getNumBins() {
			return values.size();
		}

		@Override
		public E getValue(int index) {
			return values.get(index);
		}

		@Override
		public double getInterpolatedBinIndex(E value) {
			int index = values.indexOf(value);
			Preconditions.checkState(index >= 0);
			return index;
		}
		
	}

	AbstractGMPEInterpolation(String name) {
		this.name = name;
	}
	
	public String getName() {
		return name;
	}
	
	public abstract int getNumBins();
	
	public abstract E getValue(int index);
	
	public abstract double getInterpolatedBinIndex(E value);
	
	public abstract void setGMPE_Params(ScalarIMR gmpe, ProbEqkSource source, int index);
	
//	public abstract Set<EqkRupture> getViableRuptures(Set<EqkRupture> ruptures, int index);
	
	public abstract E detectCurrentVal(ScalarIMR gmpe, Site site);
	
	public double detectInterpolatedBinIndex(ScalarIMR gmpe, Site site) {
		E curVal = detectCurrentVal(gmpe, site);
		return getInterpolatedBinIndex(curVal);
	}

	@Override
	public Iterator<E> iterator() {
		return new Iterator<E>() {
			
			private int index = 0;

			@Override
			public boolean hasNext() {
				return index < getNumBins();
			}

			@Override
			public E next() {
				return getValue(index++);
			}
		};
	}

}
