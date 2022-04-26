package scratch.kevin.nshm23;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.SubSeisMoRateReduction;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

public class TargetMFDsSumFigure {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(
//				new File("/data/kevin/markdown/inversions/fm3_1_u3ref_uniform_reproduce_ucerf3.zip"));
				new File("/data/kevin/markdown/inversions/fm3_1_u3ref_uniform_coulomb.zip"));
		
		for (double b : new double[] {0d, 0.8d, 1d}) {
			SupraSeisBValInversionTargetMFDs.Builder builder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, b);
			
			builder.sparseGR(true);
//			builder.magDepDefaultRelStdDev(M->0.1);
			builder.magDepDefaultRelStdDev(M->0.1*Math.pow(10, b*0.5*(M-6)));
			builder.applyDefModelUncertainties(true);
			builder.addSectCountUncertainties(false);
			builder.subSeisMoRateReduction(SubSeisMoRateReduction.SUB_SEIS_B_1);
			
			SupraSeisBValInversionTargetMFDs mfds = builder.build();
			for (boolean aggParents : new boolean[] {false, true}) {
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				PlotCurveCharacterstics sectChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f,
						new Color(0, 0, 0, aggParents ? 100 : 60));
				
				double leftScalar = Math.pow(10, 0)/Math.pow(10, -b*0.05);
				double rightScalar = Math.pow(10, -b*0.1)/Math.pow(10, -b*0.05);
				System.out.println("Scalars: "+leftScalar+" "+rightScalar);
				
				IncrementalMagFreqDist total = mfds.getTotalOnFaultSupraSeisMFD();
				
				List<IncrementalMagFreqDist> sectMFDs = new ArrayList<>();
				if (aggParents) {
					Map<Integer, SummedMagFreqDist> parentIDsMap = new HashMap<>();
					List<UncertainIncrMagFreqDist> subSectMFDs = mfds.getOnFaultSupraSeisNucleationMFDs();
					for (int s=0; s<subSectMFDs.size(); s++) {
						int parent = rupSet.getFaultSectionData(s).getParentSectionId();
						SummedMagFreqDist parentMFD = parentIDsMap.get(parent);
						if (parentMFD == null) {
							parentMFD = new SummedMagFreqDist(total.getMinX(), total.size(), total.getDelta());
							parentIDsMap.put(parent, parentMFD);
						}
						parentMFD.addIncrementalMagFreqDist(subSectMFDs.get(s));
					}
					sectMFDs.addAll(parentIDsMap.values());
				} else {
					sectMFDs.addAll(mfds.getOnFaultSupraSeisNucleationMFDs());
				}
				
				double minMag = Double.POSITIVE_INFINITY;
				double maxMag = 0d;
				double minY = Double.POSITIVE_INFINITY;
				for (IncrementalMagFreqDist mfd : sectMFDs) {
					DefaultXY_DataSet curLine = null;
					
					double halfDelta = 0.5*mfd.getDelta();
					
					for (Point2D pt : mfd) {
						if (pt.getY() > 0d) {
							if (curLine == null)
								curLine = new DefaultXY_DataSet();
//							curLine.set(pt.getX()-halfDelta, pt.getY());
//							curLine.set(pt.getX()+halfDelta, pt.getY());
							
							// G-R slope
							curLine.set(pt.getX()-halfDelta, pt.getY()*leftScalar);
							curLine.set(pt.getX()+halfDelta, pt.getY()*rightScalar);
							
							minMag = Math.min(minMag, pt.getX()-halfDelta);
							maxMag = Math.max(maxMag, pt.getX()+halfDelta);
							minY = Math.min(minY, pt.getY());
						} else if (curLine != null) {
							// add it
							funcs.add(curLine);
							chars.add(sectChar);
							curLine = null;
						}
					}
					if (curLine != null) {
						funcs.add(curLine);
						chars.add(sectChar);
					}
				}
				
				Range xRange = new Range(minMag, maxMag);
				minY = Math.pow(10, Math.floor(Math.log10(minY)));
				double maxY = Math.pow(10, Math.ceil(Math.log10(total.getMaxY())));
//				Range yRange = new Range(minY, maxY);
				Range yRange = new Range(1e-10, 1e-1);
				
				if (aggParents)
					funcs.get(0).setName("Parent Section Target MFDs");
				else
					funcs.get(0).setName("Subsection Target MFDs");
				
				total.setName("Total Regional Target MFD");
				funcs.add(total);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				
				PlotSpec spec = new PlotSpec(funcs, chars, "Supra-Seis MFD Construction, b="+(float)b,
						"Magnitude", "Nucleation Rate (1/yr)");
				spec.setLegendVisible(true);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(spec, false, true, xRange, yRange);
				
				String prefix = "supra_seis_mfds_b_"+(float)b;
				if (aggParents)
					prefix += "_parents";
				
				PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
			}
		}
	}

}
