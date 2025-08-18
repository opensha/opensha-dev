package scratch.kevin.prvi25;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Ints;

public class AnegadaMFDTests {
	
	private enum MFDType {
		LEN_LIMIT("L<75 km"),
		NO_LIMIT("No Limit"),
		MAG_CORNER("Mag-Corner");
		
		private String label;

		private MFDType(String label) {
			this.label = label;
		}
	}

	public static void main(String[] args) throws IOException {
		PRVI25_InvConfigFactory factory = new PRVI25_InvConfigFactory();
		LogicTreeBranch<LogicTreeNode> branch = PRVI25_LogicTree.DEFAULT_CRUSTAL_ON_FAULT.copy();
		FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, 32);
		
		File outputDir = new File("/tmp");
		
		int[] parents = {
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Anegada", "SW PROXY"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Anegada", "NE PROXY")
		};
		
		List<Table<SupraSeisBValues, NSHM23_SegmentationModels, IncrementalMagFreqDist>> typeMFDTables = new ArrayList<>();
		
		List<List<IncrementalMagFreqDist>> incrFuncLists = new ArrayList<>();
		List<List<PlotCurveCharacterstics>> charLists = new ArrayList<>();
		List<IncrementalMagFreqDist> averages = new ArrayList<>();
		MFDType[] types = MFDType.values();
		PlotLineType[] avgLineTypes = {
				PlotLineType.DOTTED,
				PlotLineType.SHORT_DASHED,
				PlotLineType.DASHED
		};
		
		CPT segCPT = GMT_CPT_Files.CATEGORICAL_TAB10_NOGRAY.instance();
		CPT bCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 1d);
		
		for (MFDType type : types) {
			List<IncrementalMagFreqDist> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			if (type == MFDType.LEN_LIMIT) {
				PRVI25_InvConfigFactory.MAX_PROXY_FAULT_RUP_LEN = 75d;
			} else {
				PRVI25_InvConfigFactory.MAX_PROXY_FAULT_RUP_LEN = 0d;
			}
			
			double avgSumWeight = 0d;
			SummedMagFreqDist averageMFD = null;
			Table<SupraSeisBValues, NSHM23_SegmentationModels, IncrementalMagFreqDist> mfdTable = HashBasedTable.create();
			typeMFDTables.add(mfdTable);
			for (SupraSeisBValues b : SupraSeisBValues.values()) {
				if (b.weight == 0d)
					continue;
				Color color = bCPT.getColor((float)b.bValue);
				double avgForBSumWeight = 0d;
				SummedMagFreqDist averageMFDforB = null;
				for (NSHM23_SegmentationModels seg : NSHM23_SegmentationModels.values()) {
					branch.setValue(b);
					double segWeight = seg.getNodeWeight(branch);
					if (segWeight == 0d)
						continue;
					if (type == MFDType.MAG_CORNER) {
						PRVI25_InvConfigFactory.PROXY_FAULT_FORCE_CLASSIC_B_1 = true;
						PRVI25_InvConfigFactory.PROXY_FAULT_MAG_CORNERS = true;
						PRVI25_InvConfigFactory.PROXY_FAULT_CLASSIC_MMAX = 7.5;
					} else {
						PRVI25_InvConfigFactory.PROXY_FAULT_FORCE_CLASSIC_B_1 = false;
						PRVI25_InvConfigFactory.PROXY_FAULT_MAG_CORNERS = false;
						PRVI25_InvConfigFactory.PROXY_FAULT_CLASSIC_MMAX = Double.NaN;
						if (seg != NSHM23_SegmentationModels.NONE)
							continue;
					}
					branch.setValue(seg);
					SummedMagFreqDist mfd = null;
					double sumWeight = 0d;
					for (NSHM23_ScalingRelationships scale : NSHM23_ScalingRelationships.values()) {
						double weight = scale.getNodeWeight(branch);
						if (weight == 0d)
							continue;
						sumWeight += weight;
						branch.setValue(scale);
						rupSet = factory.updateRuptureSetForBranch(rupSet, branch);
						factory.buildInversionConfig(rupSet, branch, 32);
						SummedMagFreqDist scaleMFD = getParentNuclMFD(parents, rupSet);
						scaleMFD.scale(weight);
						if (mfd == null)
							mfd = scaleMFD;
						else
							mfd.addIncrementalMagFreqDist(scaleMFD);
					}
					if (sumWeight != 1d)
						mfd.scale(1d/sumWeight);
					
					if (type == MFDType.MAG_CORNER) {
						mfd.setName(null);
						funcs.add(mfd);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(color.getRed(), color.getGreen(), color.getBlue(), 127)));
						mfdTable.put(b, seg, mfd);
					}
					IncrementalMagFreqDist weightedOverall = mfd.deepClone();
					weightedOverall.scale(b.weight*segWeight);
					avgSumWeight += b.weight*segWeight;
					if (averageMFD == null)
						averageMFD = new SummedMagFreqDist(mfd.getMinX(), mfd.getMaxX(), mfd.size());
					averageMFD.addIncrementalMagFreqDist(weightedOverall);
					IncrementalMagFreqDist weightedForB = mfd.deepClone();
					weightedForB.scale(segWeight);
					avgForBSumWeight += segWeight;
					if (averageMFDforB == null)
						averageMFDforB = new SummedMagFreqDist(mfd.getMinX(), mfd.getMaxX(), mfd.size());
					averageMFDforB.addIncrementalMagFreqDist(weightedForB);
				}
				
				averageMFDforB.scale(1d/avgForBSumWeight);
				averageMFDforB.setName("b="+(float)b.bValue);
				funcs.add(averageMFDforB);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2.5f, color));
				mfdTable.put(b, NSHM23_SegmentationModels.AVERAGE, averageMFDforB);
			}
			
			averageMFD.scale(1d/avgSumWeight);
			
			averageMFD.setName(type.label+" Average");
			
			incrFuncLists.add(funcs);
			charLists.add(chars);
			averages.add(averageMFD);
		}
		
		for (int t=0; t<types.length; t++) {
			MFDType type = types[t];
			
			boolean[] bySegs;
			if (type == MFDType.MAG_CORNER) {
				bySegs = new boolean[] {false, true};
			} else {
				bySegs = new boolean[] {false};
			}
			
			Table<SupraSeisBValues, NSHM23_SegmentationModels, IncrementalMagFreqDist> mfdTable = typeMFDTables.get(t);
			
			System.out.println("Plotting "+type+", table has "+mfdTable.size()+" entries");
			
			for (boolean bySeg : bySegs) {
				List<IncrementalMagFreqDist> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				if (bySeg) {
					for (NSHM23_SegmentationModels seg : NSHM23_SegmentationModels.values()) {
						if (seg == NSHM23_SegmentationModels.AVERAGE)
							continue;
						Color color = segCPT.get(seg.ordinal() % segCPT.size()).minColor;
						SummedMagFreqDist segAvg = null;
						double weightSum = 0d;
						for (SupraSeisBValues b : SupraSeisBValues.values()) {
							IncrementalMagFreqDist mfd = mfdTable.get(b, seg);
							if (mfd == null)
								continue;
							funcs.add(mfd);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(color.getRed(), color.getGreen(), color.getBlue(), 127)));
							weightSum += b.weight;
							
							mfd = mfd.deepClone();
							mfd.scale(b.weight);
							if (segAvg == null)
								segAvg = new SummedMagFreqDist(mfd.getMinX(), mfd.getMaxX(), mfd.size());
							segAvg.addIncrementalMagFreqDist(mfd);
						}
						if (segAvg == null)
							continue;
						segAvg.setName(seg.getShortName());
						segAvg.scale(1d/weightSum);
						funcs.add(segAvg);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
					}
				} else {
					for (SupraSeisBValues b : SupraSeisBValues.values()) {
						Color color = bCPT.getColor((float)b.bValue);
						for (NSHM23_SegmentationModels seg : NSHM23_SegmentationModels.values()) {
							IncrementalMagFreqDist mfd = mfdTable.get(b, seg);
							if (mfd == null)
								continue;
							if (seg == NSHM23_SegmentationModels.AVERAGE) {
								mfd.setName("b="+(float)b.bValue);
								chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
							} else {
								chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(color.getRed(), color.getGreen(), color.getBlue(), 127)));
							}
							funcs.add(mfd);
						}
					}
				}
				
				for (int i=0; i<averages.size(); i++) {
					if (i == t)
						continue;
					funcs.add(averages.get(i));
					chars.add(new PlotCurveCharacterstics(avgLineTypes[i % types.length], 3f, Color.GRAY));
				}
				funcs.add(averages.get(t));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				
				PlotSpec plot = new PlotSpec(funcs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
				plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(plot, false, true, new Range(6.5d, 8d), new Range(1e-5, 1e-1));
				
				String prefix = "anegada_mfds_"+type.name();
				if (bySeg)
					prefix += "_seg";
				
				PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, false, false);
				
				List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
				for (IncrementalMagFreqDist mfd : funcs)
					cmlFuncs.add(mfd.getCumRateDistWithOffset());
				
				plot = new PlotSpec(cmlFuncs, chars, " ", "Magnitude", "Cumulative Rate (1/yr)");
				plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				gp.drawGraphPanel(plot, false, true, new Range(6.5d, 8d), new Range(1e-5, 1e-1));
				
				PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 750, true, false, false);
			}
		}
	}
	
	private static SummedMagFreqDist getParentNuclMFD(int[] parents, FaultSystemRupSet rupSet) {
		SupraSeisBValInversionTargetMFDs mfds = rupSet.requireModule(SupraSeisBValInversionTargetMFDs.class);
		
		SummedMagFreqDist ret = null;
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			if (Ints.contains(parents, sect.getParentSectionId())) {
//			if (sect.getParentSectionId() == parent) {
				IncrementalMagFreqDist mfd = mfds.getOnFaultSupraSeisNucleationMFDs().get(sect.getSectionId());
				if (ret == null) {
					ret = new SummedMagFreqDist(mfd.getMinX(), mfd.getMaxX(), mfd.size());
				} else {
					Preconditions.checkState(ret.size() == mfd.size());
					Preconditions.checkState(ret.getMinX() == mfd.getMinX());
				}
				ret.addIncrementalMagFreqDist(mfd);
			}
		}
		
		return ret;
	}

}
