package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

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

		EvenlySpacedDouble(String name, double min, double max, int numBins) {
			super(name);
			this.min = min;
			this.max = max;
			this.numBins = numBins;
			Preconditions.checkArgument(min <= max, "Min must be <= max");
			if (min < max) {
				this.delta = (max - min) / (numBins - 1);
			} else {
				// min == max
				Preconditions.checkArgument(numBins == 1, "Num bins must be 1 if min == max");
				this.delta = 0;
			}
		}

		@Override
		public int getNumBins() {
			return numBins;
		}

		@Override
		public Double getValue(int index) {
			if (index < 0 || index > ( numBins -1 ))
				throw new IndexOutOfBoundsException("no point at index "+index);
			return min + delta*index;
		}
		
		private int getClosestIndex(double value) {
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
			double beforeDelta = value - valueBefore;
			double afterDelta = valueAfter - value;
			Preconditions.checkState(beforeDelta >= 0, "Bad beforeDelta=%s for %s with range [%s %s], binBefore=%s",
					beforeDelta, value, min, max, binBefore);
			Preconditions.checkState(afterDelta >= 0, "Bad afterDelta=%s for %s with range [%s %s], binBefore=%s",
					afterDelta, value, min, max, binBefore);
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
	
	public static abstract class LogSpacedDouble extends AbstractGMPEInterpolation<Double> {
		
		private boolean includeZero;
		
		private EvenlySpacedDouble interp;

		LogSpacedDouble(String name, double minNonZero, double max, int numBins, boolean includeZero) {
			super(name);
			Preconditions.checkState(minNonZero > 0);
			Preconditions.checkState(max > minNonZero);
			this.includeZero = includeZero;
			
			interp = new EvenlySpacedDouble(name, Math.log(minNonZero), Math.log(max), numBins) {
				
				@Override
				public void setGMPE_Params(ScalarIMR gmpe, ProbEqkSource source, int index) {
					LogSpacedDouble.this.setGMPE_Params(gmpe, source, index);
				}
				
				@Override
				public Double detectCurrentVal(ScalarIMR gmpe, Site site) {
					return LogSpacedDouble.this.detectCurrentVal(gmpe, site);
				}
			};
		}

		@Override
		public int getNumBins() {
			int bins = interp.getNumBins();
			if (includeZero)
				bins++;
			return bins;
		}

		@Override
		public Double getValue(int index) {
			if (includeZero) {
				if (index == 0)
					return 0d;
				return Math.exp(interp.getValue(index-1));
			}
			return Math.exp(interp.getValue(index));
		}

		@Override
		public double getInterpolatedBinIndex(Double value) {
			if (includeZero) {
				Preconditions.checkState(value >= 0d);
				if (value == 0d)
					return 0d;
				double logValue = Math.log(value);
				if (logValue < interp.getMin()) {
					// we're between zero and the minimum value of the log interpolator
					// do a linear interpolation between zero and the minimum value
					double index = value/Math.exp(interp.getMin());
					Preconditions.checkState(index >= 0 && index <= 1);
					return index;
				} else {
					return 1d + interp.getInterpolatedBinIndex(logValue);
				}
			}
			Preconditions.checkState(value > 0d);
			return interp.getInterpolatedBinIndex(Math.log(value));
		}
		
		public double getMin() {
			if (includeZero)
				return 0d;
			return Math.exp(interp.getMin());
		}
		
		public double getMax() {
			return Math.exp(interp.getMax());
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
