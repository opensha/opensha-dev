package scratch.kevin.nshm27;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.IntSummaryStatistics;
import java.util.List;
import java.util.concurrent.CompletableFuture;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.gui.plot.AnimatedGIFRenderer;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.RectangularDownDipSubductionRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.RupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.downDip.RectangularDownDipGrowingStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.downDip.RectangularDownDipGrowingStrategy.NeighborOverlaps;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;

import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_InterfaceFaultModels;
import net.mahdilamb.colormap.Colors;

public class DownDipRupSetBuildingTests {

	public static void main(String[] args) throws IOException {
		File baseOutputDir = new File("/home/kevin/OpenSHA/nshm26/down-dip-subsectioning");
		Preconditions.checkState(baseOutputDir.exists() || baseOutputDir.mkdir());
		
//		NSHM26_SubductionInterfaceFaultModels fm = NSHM26_SubductionInterfaceFaultModels.KERMADEC;
//		String prefix = "ker_slab2";
		
		NSHM26_InterfaceFaultModels fm = NSHM26_InterfaceFaultModels.GNMI_V1;
		String prefix = "izu_slab2";
		
		Range<Double> minSupraRange = Range.closed(20d, 40d);
		double maxSubSeisAspectRatio = 4d;
		
		File inDir = new File(baseOutputDir, prefix);
		List<? extends FaultSection> sects = fm.buildSubSects(fm);
		File outDir = new File(inDir, "rup_set_debug");
		Preconditions.checkState(outDir.exists() || outDir.mkdir());
		
		boolean doAnimation = true;
		boolean doSubSeisAnimation = false;
		boolean writeIndvFrames = true;
		boolean doMiddleRowOverlaps = false;
		boolean buildRupSet = true;
		
		FaultSubsectionCluster fullCluster = new FaultSubsectionCluster(sects);
		NeighborOverlaps neighborOverlaps = new NeighborOverlaps(fullCluster, false);
		NeighborOverlaps neighborOverlapsDD = new NeighborOverlaps(fullCluster, true);
		
		List<FaultSection> debugSects = new ArrayList<>();
		debugSects.add(sects.get(276));
		
		int numRows = sects.stream().mapToInt(s-> s.getSubSectionIndexDownDip()).max().getAsInt()+1;
		List<List<FaultSection>> rowColOrganized = new ArrayList<>(numRows);
		for (int i=0; i<numRows; i++)
			rowColOrganized.add(new ArrayList<>());
		for (FaultSection sect : sects) {
			int row = sect.getSubSectionIndexDownDip();
			int col = sect.getSubSectionIndexAlong();
			List<FaultSection> rowSects = rowColOrganized.get(row);
			while (rowSects.size() <= col)
				rowSects.add(null);
			rowSects.set(col, sect);
		}
		for (int row=0; row<numRows; row++)
			System.out.println("Row "+row+" has "+rowColOrganized.get(row).size()+" sects");
		
		// add corners
		List<FaultSection> topRow = rowColOrganized.get(0);
		List<FaultSection> bottomRow = rowColOrganized.get(numRows-1);
		debugSects.add(topRow.get(0));
		debugSects.add(topRow.get(topRow.size()-1));
		debugSects.add(bottomRow.get(0));
		debugSects.add(bottomRow.get(bottomRow.size()-1));
		
		List<FaultSection> middleRow = rowColOrganized.get(numRows/2);
		debugSects.add(middleRow.get(0));
		for (int i=1; i<5; i++)
			debugSects.add(middleRow.get((int)(middleRow.size() * i/5d)));
		debugSects.add(middleRow.get(middleRow.size()-1));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(sects);
		
		mapMaker.setFillSurfaces(true);
		mapMaker.setWriteGeoJSON(false);
		
		Color startSectColor = Colors.tab_blue;
		CPT overlapCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().reverse().rescale(0d, 1d);
		overlapCPT.setNanColor(startSectColor);
		
		Color participatingColor = overlapCPT.getMaxColor();
		Color otherColor = overlapCPT.getMinColor();
		
		RectangularDownDipGrowingStrategy growingStrat = new RectangularDownDipGrowingStrategy(minSupraRange, maxSubSeisAspectRatio);
		
		RupSetScalingRelationship scale = PRVI25_SubductionScalingRelationships.LOGA_C4p0;
		HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
		int gifWidth = 800;
		double gifFPS = 5;
		
		if (doMiddleRowOverlaps) {
			// do middle row down-dip
			File middleRowDir = new File(outDir, "middle_row_overlap_dd");
			Preconditions.checkArgument(middleRowDir.exists() || middleRowDir.mkdir());
			for (FaultSection sect : middleRow) {
				mapMaker.plotSectScalars(s-> s == sect ? Double.NaN : neighborOverlapsDD.getOverlap(sect, s),
						overlapCPT, "Fractional Overlap");
				
				mapMaker.plot(middleRowDir, "col"+frameDF.format(sect.getSubSectionIndexAlong())+"_sect"+sect.getSectionId(), " ");
				
				if (sect.getSectionId() == 361 || sect.getSectionId() == 362) {
					System.out.println("DD overlap debug from "+sect.getSectionName());
					for (FaultSection oSect : sects) {
						double overlap = neighborOverlapsDD.getOverlap(sect, oSect);
						if (overlap > 0)
							System.out.println("\t"+oSect.getSectionName()+":\t"+overlap);
					}
				}
			}
		}
		
		File animDir = new File(outDir, "animations");
		Preconditions.checkArgument(!doAnimation || animDir.exists() || animDir.mkdir());
		
		for (FaultSection debugSect : debugSects) {
			String sectPrefix = prefix+"_sect"+debugSect.getSectionId();
			System.out.println("Processing "+debugSect);
			mapMaker.plotSectScalars(s-> s == debugSect ? Double.NaN : neighborOverlaps.getOverlap(debugSect, s),
					overlapCPT, "Fractional Overlap");
			
			mapMaker.plot(outDir, sectPrefix+"_overlap", " ");
			mapMaker.plotSectScalars(s-> s == debugSect ? Double.NaN : neighborOverlapsDD.getOverlap(debugSect, s),
					overlapCPT, "Fractional Overlap");
			
			mapMaker.plot(outDir, sectPrefix+"_overlap_dd", " ");
			
			System.out.println("\tBuilding variations");
			growingStrat.setDebug(debugSect == debugSects.get(0));
			List<FaultSubsectionCluster> variations = growingStrat.getVariations(fullCluster, debugSect);
			System.out.println("\tBuilt "+variations.size()+" variations");
			Preconditions.checkState(!variations.isEmpty());
			mapMaker.clearSectScalars();
			
			if (doAnimation) {
				System.out.println("\tMaking animation");
				double minMag = Double.POSITIVE_INFINITY;
				double maxMag = 0d;
				
				AnimatedGIFRenderer gif = new AnimatedGIFRenderer(new File(animDir, sectPrefix+"_rups.gif"), gifFPS, true);
				AnimatedGIFRenderer gifSubSeis = doSubSeisAnimation ?
						new AnimatedGIFRenderer(new File(animDir, sectPrefix+"_rups_sub_seis.gif"), 0.5, true) : null;

				File fullFrameDir = new File(outDir, sectPrefix+"_rups");
				Preconditions.checkState(!writeIndvFrames || fullFrameDir.exists() || fullFrameDir.mkdir());
				
				CompletableFuture<Void> gifWriteFuture = null;
				CompletableFuture<Void> subSeisFinalizeFuture = null;
				boolean hasFullWidth = false;
				for (int i=0; i<variations.size(); i++) {
					FaultSubsectionCluster cluster = variations.get(i);
					List<Color> sectColors = new ArrayList<>(sects.size());
					for (FaultSection sect : sects) {
						if (sect == debugSect)
							sectColors.add(startSectColor);
						else if (cluster.contains(sect))
							sectColors.add(participatingColor);
						else
							sectColors.add(otherColor);
					}
					mapMaker.plotSectColors(sectColors);
					
					// get it organized
					double area = 0d;
					for (FaultSection sect : cluster.subSects)
						area += sect.getArea(true);
					
					double mag = scale.getMag(area, Double.NaN, Double.NaN, Double.NaN, cluster.startSect.getAveRake());
					minMag = Math.min(minMag, mag);
					maxMag = Math.max(maxMag, mag);
					
					String title = "Rupture "+i+"/"+variations.size()+", M"+twoDF.format(mag);
					
					PlotSpec plot = mapMaker.buildPlot(title);
					
					gp.drawGraphPanel(plot, false, false, mapMaker.getXRange(), mapMaker.getYRange());
					
					double tick = mapMaker.getAxisTick();
					PlotUtils.setXTick(gp, tick);
					PlotUtils.setYTick(gp, tick);
					
					int height = PlotUtils.calcHeight(gp, gifWidth, true);
					
					gp.getChartPanel().setSize(gifWidth, height);
					BufferedImage img = gp.getBufferedImage(gifWidth, height);
					
					if (writeIndvFrames) {
						File frameFile = new File(fullFrameDir, "rupture_"+frameDF.format(i)+".png");
						gp.saveAsPNG(frameFile.getAbsolutePath(),gifWidth, height);
					}
					
					if (gifWriteFuture != null)
						gifWriteFuture.join();
					if (hasFullWidth || !doSubSeisAnimation) {
						gifWriteFuture = CompletableFuture.runAsync(new WriteFrameRunnable(img, gif));
					} else {
						gifWriteFuture = CompletableFuture.allOf(CompletableFuture.runAsync(new WriteFrameRunnable(img, gif)),
								CompletableFuture.runAsync(new WriteFrameRunnable(img, gifSubSeis)));
					}
					
					if (!hasFullWidth && doSubSeisAnimation) {
						// see if this was full-width
						IntSummaryStatistics rowStats = cluster.subSects.stream().mapToInt(s->s.getSubSectionIndexDownDip()).summaryStatistics();
						if (rowStats.getMin() == 0 && rowStats.getMax() == rowColOrganized.size()-1) {
							// have a full width rupture
							hasFullWidth = true;
							subSeisFinalizeFuture = gifWriteFuture.thenRun(new Runnable() {
								
								@Override
								public void run() {
									try {
										gifSubSeis.finalizeAnimation();
									} catch (IOException e) {
										throw ExceptionUtils.asRuntimeException(e);
									}
								}
							});
						}
					}
				}
				gifWriteFuture.join();
				gif.finalizeAnimation();
				if (doSubSeisAnimation) {
					if (subSeisFinalizeFuture == null)
						gifSubSeis.finalizeAnimation();
					else
						subSeisFinalizeFuture.join();	
				}
				System.out.println("\tDone with animation");
			}
		}
		
		if (buildRupSet) {
			RupSetConfig rsConfig = new RectangularDownDipSubductionRupSetConfig(sects, scale);
			FaultSystemRupSet rupSet = rsConfig.build(FaultSysTools.defaultNumThreads());
			rupSet.write(new File(outDir, "rupture_set.zip"));
			File reportDir = new File(outDir, "rupture_set_report");
			
			ReportPageGen report = new ReportPageGen(rupSet, null, sects.get(0).getParentSectionName()+" Test Rupture Set",
					reportDir, ReportPageGen.getDefaultRupSetPlots(PlotLevel.DEFAULT));
			report.setReplot(true);
			report.generatePage();
		}
	}
	
	private static class WriteFrameRunnable implements Runnable {
		
		private BufferedImage img;
		private AnimatedGIFRenderer gif;

		public WriteFrameRunnable(BufferedImage img, AnimatedGIFRenderer gif) {
			this.img = img;
			this.gif = gif;
		}

		@Override
		public void run() {
			try {
				gif.writeFrame(img);
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
	}
	
	private static final DecimalFormat frameDF = new DecimalFormat("0000");
	private static final DecimalFormat twoDF = new DecimalFormat("0.00");

}
