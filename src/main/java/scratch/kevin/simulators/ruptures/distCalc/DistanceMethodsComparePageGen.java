package scratch.kevin.simulators.ruptures.distCalc;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.dom4j.DocumentException;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.distCalc.SimEventCumDistFuncSurface;
import org.opensha.sha.simulators.distCalc.SimRuptureDistCalcUtils;
import org.opensha.sha.simulators.distCalc.SimRuptureDistCalcUtils.DistanceType;
import org.opensha.sha.simulators.distCalc.SimRuptureDistCalcUtils.LocationElementDistanceCache;
import org.opensha.sha.simulators.distCalc.SimRuptureDistCalcUtils.Scalar;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class DistanceMethodsComparePageGen {
	
	private static final String SECT_BASED_NAME = "U3 Area Mapped";
	private static final double MAX_DIST = 200d;

	public static void main(String[] args) throws IOException, DocumentException {
		File mainOutputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_4983.instance();
		
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		
		double skipYears = 5000d;
		double minMag = 6.5d;
		
		// USC site
		Location siteLoc = new Location(34.0192, -118.286);
		Site site = new Site(siteLoc);
		site.addParameterList(gmpe.getSiteParams());
		LocationElementDistanceCache siteDistCache = SimRuptureDistCalcUtils.buildSiteLocDistCache(siteLoc);
		
		List<Scalar> scalars = new ArrayList<>();
		List<double[]> thresholds = new ArrayList<>();
		
		scalars.add(null); // section based
		thresholds.add(new double[] { 0d, 0.05, 0.2, 0.5 });
		
		scalars.add(Scalar.MOMENT);
		thresholds.add(new double[] { 0d, 0.01, 0.05, 0.1 });
		
//		scalars.add(Scalar.AREA);
//		thresholds.add(new double[] { 0d, 0.01, 0.05, 0.1 });
		
		File catalogOutputDir = new File(mainOutputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File outputDir = new File(catalogOutputDir, "dist_method_comparisons");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<RSQSimEvent> events = catalog.loader().skipYears(skipYears).minMag(minMag).load();
		
		System.out.println("Found "+events.size()+" events with M>"+minMag);
		
		List<String> lines = new ArrayList<>();
		lines.add("# "+catalog.getName()+" Distance Calc Comparisons");
		lines.add("");
		lines.add("This page compares various methods (each with multiple threshold values) for computing GMPE distance "
				+ "metrics from an arbitrarily complex rupture. Distances are calculated to a site at USC (*"
				+(float)siteLoc.getLatitude()+", "+(float)siteLoc.getLongitude()+"*).");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		
		for (DistanceType type : DistanceType.values()) {
			lines.add("## "+type.htmlName+" Comparisons");
			lines.add(topLink); lines.add("");
			System.out.println("Doing "+type.name());
			
			List<List<double[]>> scalarThresholdValues = new ArrayList<>();
			
			// calculate for each rupture, type, and threshold
			for (int i=0; i<scalars.size(); i++) {
				Scalar scalar = scalars.get(i);
				List<double[]> thresholdValues = new ArrayList<>();
				scalarThresholdValues.add(thresholdValues);
				for (double threshold : thresholds.get(i)) {
					if (scalar == null)
						mapper.setMinFractForInclusion(threshold);
					double[] dists = new double[events.size()];
					thresholdValues.add(dists);
					System.out.println("Calculating for "+(scalar == null ? SECT_BASED_NAME : scalar.displayName)
							+", "+(float)threshold);
					for (int j=0; j<dists.length; j++) {
						RSQSimEvent event = events.get(j);
						if (scalar == null) {
							// sub section
							RSQSimSubSectEqkRupture mapped = RSQSimUtils.buildSubSectBasedRupture(mapper, event);
							switch (type) {
							case R_JB:
								dists[j] = mapped.getRuptureSurface().getDistanceJB(siteLoc);
								break;
							case R_RUP:
								dists[j] = mapped.getRuptureSurface().getDistanceRup(siteLoc);
								break;
							case R_SEIS:
								dists[j] = mapped.getRuptureSurface().getDistanceSeis(siteLoc);
								break;

							default:
								throw new IllegalStateException("Not yet implemented");
							}
						} else {
							DiscretizedFunc distScalarFunc = SimRuptureDistCalcUtils.calcDistScalarFunc(
									event, siteLoc, siteDistCache, type, scalar);
							if (distScalarFunc == null)
								continue;
							double maxVal = distScalarFunc.getMaxY();
							double targetVal = threshold*maxVal;
							dists[j] = Double.NaN;
							for (Point2D pt : distScalarFunc) {
								if (pt.getY() >= targetVal) {
									dists[j] = pt.getX();
									break;
								}
							}
						}
					}
				}
			}
			
			// build plot tables
			for (int s1=0; s1<scalars.size(); s1++) {
				Scalar scalar1 = scalars.get(s1);
				String name1 = scalar1 == null ? SECT_BASED_NAME : scalar1.htmlName;
				double[] thresh1 = thresholds.get(s1);
				for (int s2=s1+1; s2<scalars.size(); s2++) {
					Scalar scalar2 = scalars.get(s2);
					String name2 = scalar2 == null ? SECT_BASED_NAME : scalar2.htmlName;
					double[] thresh2 = thresholds.get(s2);
					
					System.out.println("Building plot tables for "+name1+" vs "+name2);
					
					TableBuilder scatterTable = MarkdownUtils.tableBuilder();
					TableBuilder histTable = MarkdownUtils.tableBuilder();
					
					scatterTable.initNewLine();
					scatterTable.addColumn("");
					histTable.initNewLine();
					histTable.addColumn("");
					for (int i1=0; i1<thresh1.length; i1++) {
						String title = "**"+name1+", Thresh="+(float)thresh1[i1]+"**";
						scatterTable.addColumn(title);
						histTable.addColumn(title);
					}
					scatterTable.finalizeLine();
					histTable.finalizeLine();
					
					for (int i2=0; i2<thresh2.length; i2++) {
						String title = "**"+name2+", Thresh="+(float)thresh2[i2]+"**";
						scatterTable.initNewLine();
						scatterTable.addColumn(title);
						histTable.initNewLine();
						histTable.addColumn(title);
						
						for (int i1=0; i1<thresh1.length; i1++) {
							double[] vals1 = scalarThresholdValues.get(s1).get(i1);
							double[] vals2 = scalarThresholdValues.get(s2).get(i2);
							File scatterPlot = plotScatter(resourcesDir, scalar1, thresh1[i1], vals1,
									scalar2, thresh2[i2], vals2, type);
							File histPlot = plotHist(resourcesDir, scalar1, thresh1[i1], vals1,
									scalar2, thresh2[i2], vals2, type);
							
							scatterTable.addColumn("![Scatter](resources/"+scatterPlot.getName()+")");
							histTable.addColumn("![Hist](resources/"+histPlot.getName()+")");
						}
						scatterTable.finalizeLine();
						histTable.finalizeLine();
					}
					
					lines.add("### "+type.htmlName+", "+name1+" vs "+name2+" Scatters");
					lines.add(topLink); lines.add("");
					lines.addAll(scatterTable.build());
					lines.add("");
					
					lines.add("### "+type.htmlName+", "+name1+" vs "+name2+" Histograms");
					lines.add(topLink); lines.add("");
					lines.addAll(histTable.build());
					lines.add("");
				}
			}
		}
		
		lines.add("## "+gmpe.getShortName()+" 1s SA Comparisons");
		lines.add(topLink); lines.add("");
		System.out.println("Doing "+gmpe.getShortName());
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), 1d);
		
		List<List<double[]>> scalarThresholdValues = new ArrayList<>();
		
		// calculate for each rupture, type, and threshold
		for (int i=0; i<scalars.size(); i++) {
			Scalar scalar = scalars.get(i);
			List<double[]> thresholdValues = new ArrayList<>();
			scalarThresholdValues.add(thresholdValues);
			for (double threshold : thresholds.get(i)) {
				if (scalar == null)
					mapper.setMinFractForInclusion(threshold);
				double[] vals = new double[events.size()];
				thresholdValues.add(vals);
				System.out.println("Calculating for "+(scalar == null ? SECT_BASED_NAME : scalar.displayName)
						+", "+(float)threshold);
				for (int j=0; j<vals.length; j++) {
					RSQSimEvent event = events.get(j);
					RSQSimSubSectEqkRupture subSectRup = RSQSimUtils.buildSubSectBasedRupture(mapper, event);;
					RuptureSurface surf;
					if (scalar == null) {
						// sub sections
						surf = subSectRup.getRuptureSurface();
					} else {
						double xThreshold;
						if (threshold >= 0.5)
							xThreshold = threshold;
						else
							xThreshold = Math.min(0.5d, 2*threshold);
						surf = new SimEventCumDistFuncSurface(
								event, scalar, threshold, 0d, xThreshold);
					}
					EqkRupture rup = new EqkRupture(event.getMagnitude(), subSectRup.getAveRake(), surf, null);
					gmpe.setAll(rup, site, gmpe.getIntensityMeasure());
					vals[j] = gmpe.getMean();
				}
			}
		}
		
		// build plot tables
		for (int s1=0; s1<scalars.size(); s1++) {
			Scalar scalar1 = scalars.get(s1);
			String name1 = scalar1 == null ? SECT_BASED_NAME : scalar1.htmlName;
			double[] thresh1 = thresholds.get(s1);
			for (int s2=s1+1; s2<scalars.size(); s2++) {
				Scalar scalar2 = scalars.get(s2);
				String name2 = scalar2 == null ? SECT_BASED_NAME : scalar2.htmlName;
				double[] thresh2 = thresholds.get(s2);
				
				System.out.println("Building plot tables for "+name1+" vs "+name2);
				
				TableBuilder scatterTable = MarkdownUtils.tableBuilder();
				TableBuilder histTable = MarkdownUtils.tableBuilder();
				
				scatterTable.initNewLine();
				scatterTable.addColumn("");
				histTable.initNewLine();
				histTable.addColumn("");
				for (int i1=0; i1<thresh1.length; i1++) {
					String title = "**"+name1+", Thresh="+(float)thresh1[i1]+"**";
					scatterTable.addColumn(title);
					histTable.addColumn(title);
				}
				scatterTable.finalizeLine();
				histTable.finalizeLine();
				
				for (int i2=0; i2<thresh2.length; i2++) {
					String title = "**"+name2+", Thresh="+(float)thresh2[i2]+"**";
					scatterTable.initNewLine();
					scatterTable.addColumn(title);
					histTable.initNewLine();
					histTable.addColumn(title);
					
					for (int i1=0; i1<thresh1.length; i1++) {
						double[] vals1 = scalarThresholdValues.get(s1).get(i1);
						double[] vals2 = scalarThresholdValues.get(s2).get(i2);
						File scatterPlot = plotGMPEScatter(resourcesDir, scalar1, thresh1[i1], vals1,
								scalar2, thresh2[i2], vals2, gmpe);
						File histPlot = plotGMPEHist(resourcesDir, scalar1, thresh1[i1], vals1,
								scalar2, thresh2[i2], vals2, gmpe);
						
						scatterTable.addColumn("![Scatter](resources/"+scatterPlot.getName()+")");
						histTable.addColumn("![Hist](resources/"+histPlot.getName()+")");
					}
					scatterTable.finalizeLine();
					histTable.finalizeLine();
				}
				
				lines.add("### "+gmpe.getShortName()+", "+name1+" vs "+name2+" Scatters");
				lines.add(topLink); lines.add("");
				lines.addAll(scatterTable.build());
				lines.add("");
				
				lines.add("### "+gmpe.getShortName()+", "+name1+" vs "+name2+" Histograms");
				lines.add(topLink); lines.add("");
				lines.addAll(histTable.build());
				lines.add("");
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);

		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
	}
	
	private static File plotScatter(File resourcesDir, Scalar scalar1, double threshold1, double[] vals1,
			Scalar scalar2, double threshold2, double[] vals2, DistanceType type) throws IOException {
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		Preconditions.checkState(vals1.length == vals2.length, "length mismatch: %s != %s", vals1.length, vals2.length);
		
		for (int i=0; i<vals1.length; i++) {
			double v1 = vals1[i];
			double v2 = vals2[i];
			if (v1 > MAX_DIST && v2 > MAX_DIST)
				continue;
			scatter.set(v1, v2);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		
		Range range = new Range(0d, 300d);
		
		XY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(range.getLowerBound(), range.getLowerBound());
		oneToOne.set(range.getUpperBound(), range.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.GRAY));
		
		String name1 = scalar1 == null ? SECT_BASED_NAME : "Cumulative "+scalar1.displayName;
		String name2 = scalar2 == null ? SECT_BASED_NAME : "Cumulative "+scalar2.displayName;
		name1 += ", "+(float)threshold1;
		name2 += ", "+(float)threshold2;
		
		String prefix = type.name()+"_scatter";
		if (scalar1 == null)
			prefix += "_u3_mapped";
		else
			prefix += "_"+scalar1.name();
		prefix += "_"+(float)threshold1;
		if (scalar2 == null)
			prefix += "_u3_mapped";
		else
			prefix += "_"+scalar2.name();
		prefix += "_"+(float)threshold2;
		
		PlotSpec spec = new PlotSpec(funcs, chars, name1+" vs "+name2+" Scatter",
				name1+", "+type.displayName, name2+", "+type.displayName);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, false, false, range, range);
		
		File file = new File(resourcesDir, prefix+".png");
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath());
		
		return file;
	}
	
	private static File plotHist(File resourcesDir, Scalar scalar1, double threshold1, double[] vals1,
			Scalar scalar2, double threshold2, double[] vals2, DistanceType type) throws IOException {
		HistogramFunction hist = new HistogramFunction(-50d, 51, 2d);
		Range range = new Range(-50, 50);
		
		for (int i=0; i<vals1.length; i++) {
			double v1 = vals1[i];
			double v2 = vals2[i];
			if (v1 > MAX_DIST && v2 > MAX_DIST)
				continue;
			double diff = v1 - v2;
			hist.add(hist.getClosestXIndex(diff), 1d);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		String name1 = scalar1 == null ? SECT_BASED_NAME : "Cumulative "+scalar1.displayName;
		String name2 = scalar2 == null ? SECT_BASED_NAME : "Cumulative "+scalar2.displayName;
		name1 += ", "+(float)threshold1;
		name2 += ", "+(float)threshold2;
		
		String prefix = type.name()+"_hist";
		if (scalar1 == null)
			prefix += "_u3_mapped";
		else
			prefix += "_"+scalar1.name();
		prefix += "_"+(float)threshold1;
		if (scalar2 == null)
			prefix += "_u3_mapped";
		else
			prefix += "_"+scalar2.name();
		prefix += "_"+(float)threshold2;
		
		PlotSpec spec = new PlotSpec(funcs, chars, name1+" - "+name2, type.displayName+" Difference", "Count");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, false, false, range, null);
		
		File file = new File(resourcesDir, prefix+".png");
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath());
		
		return file;
	}
	
	private static File plotGMPEScatter(File resourcesDir, Scalar scalar1, double threshold1, double[] vals1,
			Scalar scalar2, double threshold2, double[] vals2, ScalarIMR gmpe) throws IOException {
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		Preconditions.checkState(vals1.length == vals2.length, "length mismatch: %s != %s", vals1.length, vals2.length);
		
		for (int i=0; i<vals1.length; i++) {
			double v1 = vals1[i];
			double v2 = vals2[i];
			if (v1 > MAX_DIST && v2 > MAX_DIST)
				continue;
			scatter.set(v1, v2);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		
		Range range = new Range(0d, 300d);
		
		XY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(range.getLowerBound(), range.getLowerBound());
		oneToOne.set(range.getUpperBound(), range.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.GRAY));
		
		String name1 = scalar1 == null ? SECT_BASED_NAME : "Cumulative "+scalar1.displayName;
		String name2 = scalar2 == null ? SECT_BASED_NAME : "Cumulative "+scalar2.displayName;
		name1 += ", "+(float)threshold1;
		name2 += ", "+(float)threshold2;
		
		String prefix = gmpe.getShortName()+"_scatter";
		if (scalar1 == null)
			prefix += "_u3_mapped";
		else
			prefix += "_"+scalar1.name();
		prefix += "_"+(float)threshold1;
		if (scalar2 == null)
			prefix += "_u3_mapped";
		else
			prefix += "_"+scalar2.name();
		prefix += "_"+(float)threshold2;
		
		PlotSpec spec = new PlotSpec(funcs, chars, name1+" vs "+name2+" Scatter",
				name1+", Ln "+gmpe.getShortName()+" Mean", name2+", Ln "+gmpe.getShortName()+" Mean");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, false, false, range, range);
		
		File file = new File(resourcesDir, prefix+".png");
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath());
		
		return file;
	}
	
	private static File plotGMPEHist(File resourcesDir, Scalar scalar1, double threshold1, double[] vals1,
			Scalar scalar2, double threshold2, double[] vals2, ScalarIMR gmpe) throws IOException {
		HistogramFunction hist = new HistogramFunction(-50d, 51, 2d);
		Range range = new Range(-50, 50);
		
		for (int i=0; i<vals1.length; i++) {
			double v1 = vals1[i];
			double v2 = vals2[i];
			if (v1 > MAX_DIST && v2 > MAX_DIST)
				continue;
			double diff = v1 - v2;
			hist.add(hist.getClosestXIndex(diff), 1d);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		String name1 = scalar1 == null ? SECT_BASED_NAME : "Cumulative "+scalar1.displayName;
		String name2 = scalar2 == null ? SECT_BASED_NAME : "Cumulative "+scalar2.displayName;
		name1 += ", "+(float)threshold1;
		name2 += ", "+(float)threshold2;
		
		String prefix = gmpe.getShortName()+"_hist";
		if (scalar1 == null)
			prefix += "_u3_mapped";
		else
			prefix += "_"+scalar1.name();
		prefix += "_"+(float)threshold1;
		if (scalar2 == null)
			prefix += "_u3_mapped";
		else
			prefix += "_"+scalar2.name();
		prefix += "_"+(float)threshold2;
		
		PlotSpec spec = new PlotSpec(funcs, chars, name1+" - "+name2,
				gmpe.getShortName()+" Ln Difference", "Count");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, false, false, range, null);
		
		File file = new File(resourcesDir, prefix+".png");
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath());
		
		return file;
	}

}
