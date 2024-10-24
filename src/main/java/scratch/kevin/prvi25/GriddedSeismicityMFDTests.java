package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.Exact;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.magdist.TaperedGR_MagFreqDist;

import com.google.common.base.Preconditions;

public class GriddedSeismicityMFDTests {
	
	private static final double FAKE_MFD_TARGET_TOTAL = 2d;
	
	private enum InputModels {
		SIMPLE_GR {
			@Override
			public EvenlyDiscretizedFunc getInputCml(double mMin, int sizeToUse, int maxSizeToBuild) {
				GutenbergRichterMagFreqDist inputCml = new GutenbergRichterMagFreqDist(mMin, maxSizeToBuild, 0.1);
				inputCml.setAllButTotMoRate(inputCml.getMinX(), inputCml.getMaxX(), 1d, 1d);
				inputCml.scale(FAKE_MFD_TARGET_TOTAL/inputCml.getY(0));
				return inputCml;
			}
		},
		BEND_UP {
			@Override
			public EvenlyDiscretizedFunc getInputCml(double mMin, int sizeToUse, int maxSizeToBuild) {
				GutenbergRichterMagFreqDist gr1 = new GutenbergRichterMagFreqDist(mMin + 0.05, maxSizeToBuild+1, 0.1);
				gr1.setAllButTotMoRate(gr1.getMinX(), gr1.getMaxX(), 9d, 1d);
				GutenbergRichterMagFreqDist gr2 = new GutenbergRichterMagFreqDist(gr1.getMinX(), gr1.size(), gr1.getDelta());
				gr2.setAllButTotMoRate(gr2.getMinX(), gr2.getMaxX(), 0.3d, 0d);
				SummedMagFreqDist summedIncr = new SummedMagFreqDist(gr1.getMinX(), gr1.size(), gr1.getDelta());
				summedIncr.addIncrementalMagFreqDist(gr1);
				summedIncr.addIncrementalMagFreqDist(gr2);
				EvenlyDiscretizedFunc inputCml = summedIncr.getCumRateDistWithOffset();
				inputCml.scale(FAKE_MFD_TARGET_TOTAL/inputCml.getY(0));
				return inputCml;
			}
		},
		BEND_DOWN {
			@Override
			public EvenlyDiscretizedFunc getInputCml(double mMin, int sizeToUse, int maxSizeToBuild) {
//				TaperedGR_MagFreqDist taper = new TaperedGR_MagFreqDist(mMin, sizeToUse+1, 0.1);
				TaperedGR_MagFreqDist taper = new TaperedGR_MagFreqDist(mMin, maxSizeToBuild, 0.1);
				taper.setAllButTotMoRate(taper.getMinX(), 7.5d, 1d, 1d);
				
//				EvenlyDiscretizedFunc inputCml = new EvenlyDiscretizedFunc(taper.getMinX(), taper.size(), taper.getDelta());
//				for (int i=0; i<taper.size(); i++)
//					inputCml.set(i, taper.getY(i));
				
				// the taper actually bends up slightly, fix that
				GutenbergRichterMagFreqDist grForMax = null;
				int taperIndex = -1;
				for (int i=0; i<taper.size(); i++) {
					GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(taper.getMinX(), taper.size(), taper.getDelta());
					gr.setAllButTotMoRate(gr.getMinX(), gr.getMaxX(), 1d, 1d);
					gr.scaleToIncrRate(taper.getX(i), taper.getY(i));
					Preconditions.checkState((float)gr.getY(i) == (float)taper.getY(i), "%s != %s, i=%s",
							(float)gr.getY(i), (float)taper.getY(i), i);
					if (grForMax == null || gr.getY(0) > grForMax.getY(0)) {
						grForMax = gr;
						taperIndex = i;
					}
				}
				System.out.println("taperIndex="+taperIndex+" at M"+(float)taper.getX(taperIndex));
				EvenlyDiscretizedFunc inputCml = new EvenlyDiscretizedFunc(taper.getMinX(), taper.size(), taper.getDelta());
				for (int i=0; i<taperIndex; i++)
					inputCml.set(i, grForMax.getY(i));
				for (int i=taperIndex; i<taper.size(); i++)
					inputCml.set(i, taper.getY(i));
				
				inputCml.scale(FAKE_MFD_TARGET_TOTAL/inputCml.getY(0));
				return inputCml;
			}
		},
		CRUSTAL_UPPER {
			@Override
			public EvenlyDiscretizedFunc getInputCml(double mMin, int sizeToUse, int maxSizeToBuild) {
				try {
					EvenlyDiscretizedFunc data = loadData("/data/erf/prvi25/seismicity/rates/2024_10_23/CRUSTAL.csv", 0.975, mMin);
//					System.out.println(data);
					return data;
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
		},
		CRUSTAL_LOWER {
			@Override
			public EvenlyDiscretizedFunc getInputCml(double mMin, int sizeToUse, int maxSizeToBuild) {
				try {
					EvenlyDiscretizedFunc data = loadData("/data/erf/prvi25/seismicity/rates/2024_10_23/CRUSTAL.csv", 0.025, mMin);
//					System.out.println(data);
					return data;
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
		};
		
		public abstract EvenlyDiscretizedFunc getInputCml(double mMin, int sizeToUse, int maxSizeToBuild);
	}
	
	private static EvenlyDiscretizedFunc loadData(String resourceName, double quantile, double mMin) throws IOException {
		CSVFile<String> csv = CSVFile.readStream(SeismicityRateFileLoader.class.getResourceAsStream(resourceName), false);
		List<Exact> branches = SeismicityRateFileLoader.loadExactBranches(csv);
		for (Exact exact : branches) {
			if ((float)exact.quantile == (float)quantile) {
				EvenlyDiscretizedFunc fullDist = exact.cumulativeDist;
				int index = fullDist.getClosestXIndex(mMin);
				if (index > 0) {
					int num = fullDist.size() - index;
					EvenlyDiscretizedFunc trimmed = new EvenlyDiscretizedFunc(
							fullDist.getX(index), num, fullDist.getDelta());
					for (int i=0; i<trimmed.size(); i++)
						trimmed.set(i, fullDist.getY(i+index));
					Preconditions.checkState((float)trimmed.getMaxX() == (float)fullDist.getMaxX());
					return trimmed;
				}
				return fullDist;
			}
		}
		throw new IllegalStateException("Didn't find quantile: %s");
	}

	public static void main(String[] args) throws IOException {
		// treat this as a cumulative
		int sizeToBuild = 91;
		boolean rescale = true;
		
		double mMin = 5d;
		double[] mMaxs = {7.9, 7.6, 7.3};
		
		EvenlyDiscretizedFunc refCml = new EvenlyDiscretizedFunc(mMin, sizeToBuild, 0.1);
		int sizeToUse = refCml.getClosestXIndex(StatUtils.max(mMaxs));
		System.out.println("SizeToUse="+sizeToUse+"; 'infinite' Mmax="+refCml.getMaxX());
		
		InputModels[] inputs = InputModels.values();
		
		for (InputModels input : inputs) {
			EvenlyDiscretizedFunc inputCml = input.getInputCml(5d, sizeToUse, sizeToBuild);
			
			System.out.println("Total rate: "+(float)inputCml.getY(0));

			//				for (boolean cmlLop : new boolean[] {false,true}) {
			for (boolean cmlLop : new boolean[] {true}) {
				
				System.out.println("Plotting "+input+", cmlLop="+cmlLop);

				List<EvenlyDiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();

				CPT cpt = GMT_CPT_Files.CATEGORICAL_TAB10.instance();

				inputCml.setName("Original (Input) Cumulative MFD");
				funcs.add(inputCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 7f, cpt.get(chars.size()).minColor));
				
				for (double mMax : mMaxs) {
					int size = refCml.getClosestXIndex(mMax);
					IncrementalMagFreqDist incr = new IncrementalMagFreqDist(
							inputCml.getMinX()+0.5*inputCml.getDelta(), size, inputCml.getDelta());
					
					if (cmlLop) {
						// convert to incremental
						for (int i=0; i<incr.size(); i++) {
							double binStart = inputCml.getY(i);
							double binEnd = i<inputCml.size()-1 ? inputCml.getY(i+1) : 0;
							incr.set(i, binStart - binEnd);
						}
					} else {
						// build full incremental
						IncrementalMagFreqDist incrGR = new IncrementalMagFreqDist(
								inputCml.getMinX()+0.5*inputCml.getDelta(), inputCml.size()-1, inputCml.getDelta());
						for (int i=0; i<incrGR.size(); i++) {
							double binStart = inputCml.getY(i);
							double binEnd = i<inputCml.size()-1 ? inputCml.getY(i+1) : 0;
							incrGR.set(i, binStart - binEnd);
						}

						for (int i=0; i<incr.size(); i++)
							incr.set(i, incrGR.getY(i));
					}
					
					if (rescale) {
						double factor = inputCml.getY(0)/incr.calcSumOfY_Vals();
						System.out.println("Scaling by "+(float)factor+" for Mmax="+(float)mMax);
						incr.scale(factor);
					}
					
					EvenlyDiscretizedFunc cml = incr.getCumRateDistWithOffset();
					
					incr.setName("Converted to Incremental (Mmax="+(float)mMax+")");
					funcs.add(incr);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 7f, cpt.get(chars.size()).minColor));

					cml.setName("Back to Cumulative (Mmax="+(float)mMax+")");
					funcs.add(cml);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, cpt.get(chars.size()).minColor));
				}

				//					String title = cmlLop ? "Cumulative Lopping" : "Incremental Lopping";
				String title = " ";
				String prefix = "mfd_lop_tests";
				if (!cmlLop)
					prefix += "_incr_lop";
				prefix += "_"+input.name();

				PlotSpec spec = new PlotSpec(funcs, chars, title, "Magnitude", "Rate");
				spec.setLegendInset(true);

				HeadlessGraphPanel gp = PlotUtils.initHeadless();

				gp.drawGraphPanel(spec, false, true, new Range(5d, 8.1d), new Range(1e-4, 10));

				PlotUtils.writePlots(new File("/tmp"), prefix, gp, 1000, 900, true, false, false);
			}
		}
	}

}
