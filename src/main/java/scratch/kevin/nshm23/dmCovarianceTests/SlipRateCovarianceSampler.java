package scratch.kevin.nshm23.dmCovarianceTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SlipRatePlots;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

public class SlipRateCovarianceSampler {
	
	private SectionCovarianceSampler sampler;
	
	private List<double[]> prevSamples;
	private List<List<FaultSection>> prevSampled;
	
	public static MeanPreservationAdjustment MEAN_ADJUSTMENT_DEFAULT = MeanPreservationAdjustment.ADJUST_MEAN;
	private MeanPreservationAdjustment meanAdjustment = MEAN_ADJUSTMENT_DEFAULT;
	
	public enum MeanPreservationAdjustment {
		CAP_SIGMA, // this preserves mean perfectly, but reduces variability when sigma is near or greater than mean
		ADJUST_MEAN // this does a better job of preserving variability, but results in mean bias when sigma is greater than mean
	}
	
	/*
	 * Sigma truncation: don't sample values more than these many std. devs. away from the mean to avoid extremem outliers
	 */
	public static double TRUNCATION_DEFAULT = 3d;
	private double truncation = TRUNCATION_DEFAULT;
	
	/*
	 * Setting many things to zero messes up target MFD calculations. Force each sampled slip rate to be at least this
	 * fraction of the original slip rate (still allows zero if the original slip rate was zero)
	 */
	public static double FRACT_FLOOR_DEFAULT = 1e-2;
	private double fractFloor = FRACT_FLOOR_DEFAULT;

	private List<? extends FaultSection> refSects;

	public SlipRateCovarianceSampler(SectionCovarianceSampler sampler) {
		this(sampler, sampler.getSubSects());
	}

	public SlipRateCovarianceSampler(SectionCovarianceSampler sampler, List<? extends FaultSection> refSects) {
		this.sampler = sampler;
		this.refSects = refSects;
	}
	
	public List<List<FaultSection>> buildSamples(int numSamples, RandomGenerator rng, boolean center, int interpSkip) {
		return buildSamples(sampler.sample(numSamples, rng, center, interpSkip));
	}
	
	public List<FaultSection> buildSample(RandomGenerator rng, int interpSkip) {
		return buildSamples(sampler.sample(1, rng, false, interpSkip)).get(0);
	}
	
	private class Incrementer {
		int val = 0;
		
		public synchronized void increment() {
			val++;
		}
		
		@Override
		public String toString() {
			return val+"";
		}
	}
	
	public List<List<FaultSection>> buildSamples(List<double[]> samples) {
		List<List<FaultSection>> ret = new ArrayList<>();
		
		int numSects = refSects.size();
		
		MinMaxAveTracker overallZTrack = new MinMaxAveTracker();
		int totNumSetToZero = 0;
		int totNumSetToFloor = 0;
		int totNumCapped = 0;
		int totNumTruncated = 0;
		int totNumZerosCappedAtOneSigma = 0;
		
		DecimalFormat pDF = new DecimalFormat("0.0%");
		
		boolean capSigma = meanAdjustment == MeanPreservationAdjustment.CAP_SIGMA;
		MinMaxAveTracker meanAdjRatioTrack = new MinMaxAveTracker();
		MinMaxAveTracker meanAdjDiffTrack = new MinMaxAveTracker();
		Incrementer numMeanAdjToZero = new Incrementer();
		int numMeanAdj = 0;
		
		double[] slipMeans = new double[numSects];
		double[] slipFloors = new double[numSects];
		double[] slipSigmas = new double[numSects];
		ExecutorService exec = null;
		List<Future<?>> slipAdjFutures = new ArrayList<>();
		for (int s=0; s<numSects; s++) {
			FaultSection sect = refSects.get(s);
			// TODO should this be creep reduced?
			slipMeans[s] = sect.getOrigAveSlipRate();
			slipSigmas[s] = sect.getOrigSlipRateStdDev();
			slipFloors[s] = fractFloor*slipMeans[s];
			if (meanAdjustment == MeanPreservationAdjustment.ADJUST_MEAN && slipMeans[s] > 0d) {
				// we can skip the adjustment if mean - trunc*sigma >0
				// if we're not applying truncation, pretend as though we are at 3 sigma
				double truncation = (Double.isFinite(this.truncation) && this.truncation > 0d) ? this.truncation : 3d;
				if (slipMeans[s] - truncation*slipSigmas[s] < 0d) {
					// we need to do it
					numMeanAdj++;
					
					if (exec == null)
						exec = Executors.newFixedThreadPool(Integer.min(16, FaultSysTools.defaultNumThreads()));
					
					int index = s;
					slipAdjFutures.add(exec.submit(new Runnable() {
						
						@Override
						public void run() {
							// TODO analytical solution?
							int iters = 100;
							int samplesPerIter = 50000;
							Random r = new Random(sect.getSectionId()*samples.size());
							double[] slipSamples = new double[samplesPerIter];
							double origMean = slipMeans[index];
							double curMean = origMean;
							double floor = slipFloors[index];
							double sigma = slipSigmas[index];
							for (int i=0; i<iters; i++) {
								for (int n=0; n<samplesPerIter; n++) {
									double z = r.nextGaussian();
									if (z < -truncation)
										z = -truncation;
									else if (z > truncation)
										z = truncation;
									slipSamples[n] = Math.max(floor, curMean + z*sigma);
								}
								double calcMean = StatUtils.mean(slipSamples);
								double ratio = calcMean / origMean;
								// should never go above the original slip rate
								curMean = Math.min(origMean, curMean / ratio);
							}
							double ratio = curMean/origMean;
							double diff = curMean - origMean;
							
							synchronized (SlipRateCovarianceSampler.this) {
//								System.out.println("Adjustment for "+sect.getSectionName());
//								System.out.println("\tOrig slip:\t"+(float)origMean+" +/- "+(float)sigma);
//								System.out.println("\tMod slip:\t"+(float)curMean+";\tratio="+(float)ratio+";\tdiff="+(float)diff);
								if (ratio < 1e-5) {
									// just set it to zero
									numMeanAdjToZero.increment();
									slipMeans[index] = 0d;
								} else {
									meanAdjRatioTrack.addValue(ratio);
									meanAdjDiffTrack.addValue(diff);
									slipMeans[index] = curMean;
								}
							}
						}
					}));
				}
			}
		}
		for (Future<?> slipAdjFuture : slipAdjFutures) {
			try {
				slipAdjFuture.get();
			} catch (InterruptedException | ExecutionException e) {
				exec.shutdown();
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		if (exec != null)
			exec.shutdown();
		
		if (meanAdjustment == MeanPreservationAdjustment.ADJUST_MEAN) {
			System.out.println("Adjusted "+numMeanAdj+" section mean slips for symmetry after forced positivity");
			System.out.println("\t"+numMeanAdjToZero+" adjusted to zero");
			System.out.println("\tRatio stats:\t"+meanAdjRatioTrack);
			System.out.println("\tDiff stats:\t"+meanAdjDiffTrack);
		}
		
		for (int n=0; n<samples.size(); n++) {
			double[] sample = samples.get(n);
			Preconditions.checkState(sample.length == numSects);
			List<FaultSection> sampledSects = new ArrayList<>(numSects);
			
			MinMaxAveTracker zTrack = new MinMaxAveTracker();
			int numSetToFloor = 0;
			int numSetToZero = 0;
			int numCapped = 0;
			int numTruncated = 0;
			int numZerosCappedAtOneSigma = 0;
			
			for (int s=0; s<sample.length; s++) {
				FaultSection modSect = refSects.get(s).clone();
				double origSlipRate = slipMeans[s];
				double origStdDev = slipSigmas[s];
				double z = sample[s];
				if (truncation > 0d && (z > truncation || z < -truncation)) {
					numTruncated++;
					z = Math.max(-truncation, z);
					z = Math.min(z, truncation);
				}
				if (capSigma && origSlipRate > 0d && z > 0d) {
					if (origSlipRate == 0d) {
						// for zero slip rate, cap it at +1 sigma (equivalent to what we do when sigma > mean)
						if (z > 1d) {
							numZerosCappedAtOneSigma++;
							z = 1d;
						}
					} else {
						double zeroZ = -origSlipRate/origStdDev;
						if (z > -zeroZ) {
							// this is a higher z score than we could go negative, cap it to keep averages in tact
							z = -zeroZ;
							numCapped++;
						}
					}
				}
					
				zTrack.addValue(z);
				if (origStdDev > 0d) {
					// we have a standard deviation, figure out new slip rate
					double randSlip = origSlipRate + z*origStdDev;
					if (randSlip < slipFloors[s]) {
						numSetToFloor++;
						randSlip = slipFloors[s];
						if (randSlip == 0d)
							numSetToZero++;
					}
					// TODO should this be creep reduced?
					modSect.setAveSlipRate(randSlip);
				}
				sampledSects.add(modSect);
			}
			if (n == 0 && sampler.isDebug()) {
				System.out.println("Random sample 1");
				System.out.println("\tz score stats: "+zTrack);
				if (truncation > 0d)
					System.out.println("\t"+numTruncated+" ("+pDF.format((double)numTruncated/(double)numSects)
							+") were truncated at sigma="+(float)truncation);
				System.out.println("\t"+numSetToFloor+" ("+pDF.format((double)numSetToFloor/(double)numSects)
							+") were set to floor values");
				System.out.println("\t"+numSetToZero+" ("+pDF.format((double)numSetToZero/(double)numSects)
							+") were set to zero because the sampled slip was negative and the original value was zero");
				if (capSigma) {
					System.out.println("\t"+numCapped+" ("+pDF.format((double)numCapped/(double)numSects)
					+") were capped for symmetry w.r.t zero");
					System.out.println("\t"+numZerosCappedAtOneSigma+" ("+pDF.format((double)numZerosCappedAtOneSigma/(double)numSects)
							+") were capped at +1 sigma due to zero mean");
				}
			}
			overallZTrack.addFrom(zTrack);
			totNumCapped += numCapped;
			totNumZerosCappedAtOneSigma += numZerosCappedAtOneSigma;
			totNumTruncated += numTruncated;
			totNumSetToZero += numSetToZero;
			totNumSetToFloor += numSetToFloor;
			ret.add(sampledSects);
		}
		
		double avgNumCapped = (double)totNumCapped/(double)samples.size();
		double avgNumZerosCappedAtOneSigma = (double)totNumZerosCappedAtOneSigma/(double)samples.size();
		double avgNumTruncated = (double)totNumTruncated/(double)samples.size();
		double avgNumSetToFloor = (double)totNumSetToFloor/(double)samples.size();
		double avgNumSetToZero = (double)totNumSetToZero/(double)samples.size();
		
		DecimalFormat avgDF = new DecimalFormat("0.0");
		
		System.out.println("Slip Sampling Stats");
		System.out.println("\tz score stats: "+overallZTrack);
		if (truncation > 0d)
			System.out.println("\t"+avgDF.format(avgNumTruncated)+" ("+pDF.format(avgNumTruncated/(double)numSects)
					+") were truncated at sigma="+(float)truncation);
		System.out.println("\t"+avgDF.format(avgNumSetToFloor)+" ("+pDF.format(avgNumSetToFloor/(double)numSects)
					+") were set to floor values");
		System.out.println("\t"+avgDF.format(avgNumSetToZero)+" ("+pDF.format(avgNumSetToZero/(double)numSects)
					+") were set to zero because the sampled slip was negative and the original value was zero");
		if (capSigma) {
			System.out.println("\t"+avgDF.format(avgNumCapped)+" ("+pDF.format(avgNumCapped/(double)numSects)
					+") were capped for symmetry w.r.t zero");
			System.out.println("\t"+avgDF.format(avgNumZerosCappedAtOneSigma)+" ("+pDF.format(avgNumZerosCappedAtOneSigma/(double)numSects)
			+") were capped at +1 sigma due to zero mean");
		}
		
		this.prevSamples = samples;
		this.prevSampled = ret;
		return ret;
	}
	
	public void buildDiagnosticPage(File outputDir, RupSetDeformationModel dm, int numSamplePlots, int numSectPlots) throws IOException {
		System.out.println("Building report in "+outputDir.getAbsolutePath());
		Preconditions.checkNotNull(prevSampled, "No samples have been built");
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<? extends FaultSection> subSects = sampler.getSubSects();
		GeographicMapMaker mapMaker = new GeographicMapMaker(subSects);
		mapMaker.setWritePDFs(false);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setScalarThickness(2f);
		
		CPT corrCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d);
		CPT zScoreCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-2d, 2d);
		zScoreCPT.setNanColor(Color.GRAY);
		double maxSlip = 0d;
		for (FaultSection sect : subSects)
			maxSlip = Math.max(sect.getOrigAveSlipRate(), maxSlip);
		double maxDiff;
		if (maxSlip > 30d) {
			maxSlip = 40d;
			maxDiff = 5d;
		} else if (maxSlip > 20d) {
			maxSlip = 30d;
			maxDiff = 5d;
		} else if (maxSlip > 15d) {
			maxSlip = 20d;
			maxDiff = 5d;
		} else if (maxSlip > 10d) {
			maxSlip = 15d;
			maxDiff = 3d;
		} else if (maxSlip > 5d) {
			maxSlip = 10d;
			maxDiff = 2d;
		} else {
			maxSlip = 5d;
			maxDiff = 1d;
		}
		CPT linearSlipCPT = SlipRatePlots.linearSlipCPT(maxSlip);
		CPT logSlipCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-3, 2);
		logSlipCPT.setNanColor(Color.GRAY);
		CPT logRatioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d);
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-maxDiff, maxDiff);
		CPT covCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 2d);
		covCPT.setNanColor(Color.GRAY);
		CPT covDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-0.2, 0.2);
		CPT covLogRatioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d);
		
		int numSects = subSects.size();
		int numSamples = prevSampled.size();
		double[] origSlips = new double[numSects];
		double[] origCOVs = new double[numSects];
		double[] sampledMeanSlips = new double[numSects];
		double[] sampledCOVs = new double[numSects];
		for (int s=0; s<numSects; s++) {
			origSlips[s] = subSects.get(s).getOrigAveSlipRate();
			origCOVs[s] = origSlips[s] > 0 ? subSects.get(s).getOrigSlipRateStdDev()/origSlips[s] : Double.NaN;
			double[] samples = new double[numSamples];
			for (int n=0; n<numSamples; n++)
				samples[n] = prevSampled.get(n).get(s).getOrigAveSlipRate();
			sampledMeanSlips[s] = StatUtils.mean(samples);
			sampledCOVs[s] = Math.sqrt(StatUtils.variance(samples))/sampledMeanSlips[s];
		}
		
		List<String> lines = new ArrayList<>();
		lines.add("# Slip Rate Sampling, "+dm.getShortName()+" ( "+prevSampled.size()+" samples)");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Average Slips");
		lines.add(topLink); lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.initNewLine();
		mapMaker.plotSectScalars(origSlips, linearSlipCPT, "Original Slip Rates (mm/yr)");
		mapMaker.plot(resourcesDir, "orig_slips", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/orig_slips.png)");
		mapMaker.plotSectScalars(log10(origSlips), logSlipCPT, "Log10(Original Slip Rates) (mm/yr)");
		mapMaker.plot(resourcesDir, "orig_slips_log", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/orig_slips_log.png)");
		table.finalizeLine();
		
		table.initNewLine();
		mapMaker.plotSectScalars(sampledMeanSlips, linearSlipCPT, "Sampled Mean Slip Rates (mm/yr)");
		mapMaker.plot(resourcesDir, "sampled_mean_slips", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/sampled_mean_slips.png)");
		mapMaker.plotSectScalars(log10(sampledMeanSlips), logSlipCPT, "Log10(Sampled Mean Rates) (mm/yr)");
		mapMaker.plot(resourcesDir, "sampled_mean_slips_log", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/sampled_mean_slips_log.png)");
		table.finalizeLine();
		
		table.initNewLine();
		mapMaker.plotSectScalars(diffs(origSlips, sampledMeanSlips), diffCPT, "Sampled Mean - Orig Slip Rates (mm/yr)");
		mapMaker.plot(resourcesDir, "sampled_mean_slips_diff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/sampled_mean_slips_diff.png)");
		mapMaker.plotSectScalars(logRatio(origSlips, sampledMeanSlips), logRatioCPT, "Log10(Sampled Mean/Orig Rates)");
		mapMaker.plot(resourcesDir, "sampled_mean_slips_ratio", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/sampled_mean_slips_ratio.png)");
		table.finalizeLine();
		
		table.initNewLine();
		mapMaker.plotSectScalars(origCOVs, covCPT, "Original Slip COV");
		mapMaker.plot(resourcesDir, "orig_cov", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/orig_cov.png)");
		mapMaker.plotSectScalars(sampledCOVs, covCPT, "Sampled Slip COV");
		mapMaker.plot(resourcesDir, "sampled_cov", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/sampled_cov.png)");
		table.finalizeLine();
		
		table.initNewLine();
		mapMaker.plotSectScalars(diffs(origCOVs, sampledCOVs), covDiffCPT, "Sampled COV - Orig COV");
		mapMaker.plot(resourcesDir, "sampled_cov_diff", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/sampled_cov_diff.png)");
		mapMaker.plotSectScalars(logRatio(origCOVs, sampledCOVs), covLogRatioCPT, "Log10(Sampled COV/Orig COV)");
		mapMaker.plot(resourcesDir, "sampled_cov_ratio", " ");
		table.addColumn("![Map]("+resourcesDir.getName()+"/sampled_cov_ratio.png)");
		table.finalizeLine();
		
		lines.addAll(table.build());
		
		if (numSamplePlots > 0) {
			if (numSamplePlots > numSamples)
				numSamplePlots = numSamples;
			
			lines.add("## Individual Samples");
			lines.add(topLink); lines.add("");
			table = MarkdownUtils.tableBuilder();
			
			for (int i=0; i<numSamplePlots; i++) {
				String prefix = "sample_"+i;
				
				double[] sampleZ = prevSamples.get(i);
				Preconditions.checkState(sampleZ.length == numSects);
				double[] sampleSlips = new double[numSects];
				List<FaultSection> sampled = prevSampled.get(i);
				for (int s=0; s<numSects; s++)
					sampleSlips[s] = sampled.get(s).getOrigAveSlipRate();
				
				table.initNewLine();
				
				table.addColumn(MarkdownUtils.boldCentered("Sample "+i));
				
				mapMaker.plotSectScalars(sampleSlips, linearSlipCPT, "Sample "+i+" Slip Rates (mm/yr)");
				mapMaker.plot(resourcesDir, prefix+"_slips", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_slips.png)");
				
				mapMaker.plotSectScalars(log10(sampleSlips), logSlipCPT, "Log10 (Sample "+i+" Slip Rates) (mm/yr)");
				mapMaker.plot(resourcesDir, prefix+"_slips_log", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_slips_log.png)");
				
				table.finalizeLine().initNewLine();
				
				mapMaker.plotSectScalars(sampleZ, zScoreCPT, "Sample "+i+" Slip Rate z-scores");
				mapMaker.plot(resourcesDir, prefix+"_z", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_z.png)");
				
				mapMaker.plotSectScalars(diffs(origSlips, sampleSlips), diffCPT, "Sample "+i+" - Orig Slip Rate (mm/yr)");
				mapMaker.plot(resourcesDir, prefix+"_diff", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_diff.png)");
				
				mapMaker.plotSectScalars(logRatio(origSlips, sampleSlips), logRatioCPT, "Log10(Sample "+i+" / Orig Slip Rate)");
				mapMaker.plot(resourcesDir, prefix+"_ratio", " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/"+prefix+"_ratio.png)");
				
				table.finalizeLine();
			}
			
			lines.addAll(table.build());
		}
		
		if (numSectPlots > 0) {
			Random rand = new Random((long)subSects.size()*numSamples);
			lines.add("## Section Samples");
			lines.add(topLink); lines.add("");
			table = MarkdownUtils.tableBuilder();
			
			table.initNewLine();
			for (int s=0; s<numSectPlots; s++) {
				int sectIndex = rand.nextInt(numSects);
				FaultSection sect = subSects.get(sectIndex);
				Region sectRegion = GeographicMapMaker.buildBufferedRegion(List.of(sect), 150d, false);
				Region calcRegion = GeographicMapMaker.buildBufferedRegion(List.of(sect), 200d, false);
				System.out.println("Calculating correlations for "+sect.getSectionName());
				double[] cors = new double[numSects];
				for (int i=0; i<numSects; i++) {
					if (i == sectIndex) {
						cors[i] = 1d;
					}else {
						FaultSection oSect = subSects.get(i);
						Location traceStart = oSect.getFaultTrace().first();
						Location traceEnd = oSect.getFaultTrace().last();
						Location middle = new Location(0.5d*(traceStart.getLatitude()+traceEnd.getLatitude()),
								0.5*(traceStart.getLongitude()+traceEnd.getLongitude()));
						if (calcRegion.contains(middle))
							cors[i] = sampler.getCorrelationCoefficient(sect, subSects.get(i));
					}
				}
				mapMaker.setRegion(sectRegion);
				mapMaker.setSectHighlights(List.of(sect), new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				mapMaker.plotSectScalars(cors, corrCPT, sect.getName()+" Correleation Coeffs");
				mapMaker.plot(resourcesDir, "corr_plot_"+s, " ");
				table.addColumn("![Map]("+resourcesDir.getName()+"/corr_plot_"+s+".png)");
			}
			table.finalizeLine();
			
			lines.addAll(table.wrap(3, 0).build());
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
		
		System.out.println("DONE Building report");
	}
	
	private double[] diffs(double[] orig, double[] mod) {
		double[] ret = new double[orig.length];
		for (int i=0; i<orig.length; i++)
			ret[i] = mod[i] - orig[i];
		return ret;
	}
	
	private double[] log10(double[] orig) {
		double[] ret = new double[orig.length];
		for (int i=0; i<orig.length; i++) {
			if (Double.isFinite(orig[i]) && orig[i] > 0d)
				ret[i] = Math.log10(orig[i]);
			else
				ret[i] = Double.NaN;
		}
		return ret;
	}
	
	private double[] logRatio(double[] orig, double[] mod) {
		double[] ret = new double[orig.length];
		for (int i=0; i<orig.length; i++)
			ret[i] = Math.log10(mod[i]/orig[i]);
		return ret;
	}

	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;
//		NSHM23_DeformationModels dm = NSHM23_DeformationModels.GEOLOGIC;
//		NSHM23_DeformationModels dm = NSHM23_DeformationModels.ZENG;
		NSHM23_DeformationModels dm = NSHM23_DeformationModels.POLLITZ;
//		NSHM23_DeformationModels dm = NSHM23_DeformationModels.EVANS;
//		NSHM23_DeformationModels dm = NSHM23_DeformationModels.SHEN_BIRD;
		
//		NSHM23_SingleStates state = NSHM23_SingleStates.UT;
//		NSHM23_SingleStates state = NSHM23_SingleStates.CA;
		NSHM23_SingleStates state = null;
//		int interpSkip = 2;
//		int interpSkip = 1;
		int interpSkip = 0;
		
		boolean ignoreCache = false;
		
		// disable the upper bound
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = 0d;
		
		File mainOutputDir = new File("/home/kevin/markdown/nshm23-misc/dm_slip_sampling");
		String dirPrefix = state == null ? "full" : state.getFilePrefix();
		String cacheDirPrefix = dirPrefix;
		if (dm == fm.getDefaultDeformationModel()) {
			cacheDirPrefix += "_"+dm.getFilePrefix();
			dirPrefix += "_"+dm.getFilePrefix();
		} else {
			cacheDirPrefix += "_"+fm.getDefaultDeformationModel().getFilePrefix();
			dirPrefix += "_"+dm.getFilePrefix();
		}
		if (interpSkip > 0) {
			dirPrefix += "_interpSkip"+interpSkip;
			cacheDirPrefix += "_interpSkip"+interpSkip;
		}
		
		List<? extends FaultSection> subSects = dm.build(fm);
		List<? extends FaultSection> origSubSects = subSects;
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(subSects);
		
		if (state != null) {
			List<FaultSection> filtered = new ArrayList<>();
			for (FaultSection sect : subSects)
				if (state.contains((GeoJSONFaultSection)sect))
					filtered.add(sect);
			System.out.println("Retained "+filtered.size()+"/"+subSects.size()+" sections for "+state);
			subSects = filtered;
		}
		
		double corrDist = 200d;
		double zeroDistCoeff = 0.95;
		double negativeCorrMaxDist = 30d;
//		SectionCovarianceSampler sampler = new SectionCovarianceSampler.QuickTestInvDistance(
//				subSects, distAzCalc, 100d, 0.9);
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
//				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
//		Preconditions.checkState(origSubSects.size() == sol.getRupSet().getNumSections());
//		SectionCovarianceSampler sampler = new PrecomputedConnectivityCorrelationSampler(
//				subSects, sol, distAzCalc, corrDist, zeroDistCoeff, negativeCorrMaxDist);
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		FaultSystemRupSet rupSet = factory.buildRuptureSet(NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT, FaultSysTools.defaultNumThreads());
		Preconditions.checkState(origSubSects.size() == rupSet.getNumSections());
		SectionCovarianceSampler sampler = new BvalAndSegConnectivityCorrelationSampler(
				subSects, rupSet, distAzCalc, corrDist, zeroDistCoeff, negativeCorrMaxDist, 0.5d, NSHM23_SegmentationModels.MID, true);
		
		dirPrefix += "_"+sampler.getSamplerPrefix();
		cacheDirPrefix += "_"+sampler.getSamplerPrefix();
		
		File outputDir = new File(mainOutputDir, dirPrefix);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File cachedOutputDir = dirPrefix.equals(cacheDirPrefix) ? outputDir : new File(mainOutputDir, cacheDirPrefix);
		Preconditions.checkState(cachedOutputDir.exists() || cachedOutputDir.mkdir());
		File cacheDir = new File(cachedOutputDir, "cache");
		Preconditions.checkState(cacheDir.exists() || cacheDir.mkdir());
		
		if (!ignoreCache) {
			// see if we already have it cached
			Optional<SectionCovarianceSampler> cached = sampler.loadCached(cacheDir, interpSkip);
			if (cached.isPresent()) {
				System.out.println("Using cached version!");
				sampler = cached.get();
			}
		}
		
		RandomGenerator rng = new Well19937c(12345l);
		boolean center = false;
		
		sampler.setDebug(true);
		SlipRateCovarianceSampler slipSampler = new SlipRateCovarianceSampler(sampler);
		
		int numSamples = 450;
		int numSamplePlot = 10;
		int numSectPlots = 15;
		
		slipSampler.buildSamples(numSamples, rng, center, interpSkip);

		sampler.writeCache(cacheDir);
		
		slipSampler.buildDiagnosticPage(outputDir, dm, numSamplePlot, numSectPlots);
	}

}
