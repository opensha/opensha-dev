package scratch.kevin.nshm23.timJonEpistemic;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_LogicTreeNode;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import com.google.common.base.Preconditions;

import scratch.kevin.ucerf3.TimJonEpistemicCurveCalc;

public class PostProcessCurveZips {

	public static void main(String[] args) throws IOException {
		File dir = new File("/project/scec_608/kmilner/nshm23/batch_inversions/"
				+ "2023_08_23-tim_jon_nshm23_v3_epistemic_calcs_downsampled");
		
		File outputDir = new File(dir, "results");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		LogicTree<?> tree = LogicTree.read(new File(dir, "logic_tree.json"));
		ZipFile zip = new ZipFile(new File(dir, "results_hazard_sites.zip"));
		
		List<String> siteNames = new ArrayList<>();
		List<Location> siteLocs = new ArrayList<>();
		List<double[]> siteVs30sMeasured = new ArrayList<>();
		List<double[]> siteVs30sInferred = new ArrayList<>();
		
		double[] periods = {0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75,
				1, 1.5, 2, 3, 4, 5, 7.5, 10};
		
		double[] percentiles = {
				0.1,
				1.3,
				2.5,
				5,
				9.25,
				16,
				33,
				50,
				67,
				84,
				90.75,
				95,
				97.5,
				98.7,
				99.9
		};
		
		AttenRelRef[] gmpeRefs = {
				AttenRelRef.ASK_2014,
				AttenRelRef.BSSA_2014,
				AttenRelRef.CB_2014,
				AttenRelRef.CY_2014
		};
		
		siteNames.add("Davis");
		siteLocs.add(new Location(38.3210, -121.4530));
		siteVs30sMeasured.add(new double[] {
				216,
				236,
				251,
				266,
				281,
				299,
				327
		});
		siteVs30sInferred.add(new double[] {
				179,
				212,
				238,
				266,
				296,
				333,
				394
		});
		
		siteNames.add("Berkeley");
		siteLocs.add(new Location(37.5216, -122.1527));
		siteVs30sMeasured.add(new double[] {
				514,
				562,
				598,
				633,
				670,
				713,
				779
		}); 
		siteVs30sInferred.add(new double[] {
				288,
				402,
				510,
				633,
				785,
				996,
				1393
		});
		
		for (int s=0; s<siteNames.size(); s++) {
			String siteName = siteNames.get(s);
			double[] vs30Measured = siteVs30sMeasured.get(s);
			double[] vs30Inferred = siteVs30sInferred.get(s);
			
			for (boolean inferred : new boolean[] {false,true}) {
				double[] vs30s = inferred ? vs30Inferred : vs30Measured;
				String typeParamValue = inferred ? Vs30_TypeParam.VS30_TYPE_INFERRED : Vs30_TypeParam.VS30_TYPE_MEASURED;
				
				for (double vs30 : vs30s) {
					String sitePrefix = siteName+"_"+(int)vs30+"_"+typeParamValue;
					
					for (double period : periods) {
						String entryName = sitePrefix;
						String outputPrefix = "curves_"+siteName;
						if (period == 0d) {
							entryName += "_pga";
							outputPrefix += "_pga";
						} else {
							entryName += "_sa_"+(float)period;
							outputPrefix += "_sa_"+(float)period+"s";
						}
						entryName += ".csv";
						
						System.out.println("Reading "+entryName);
						
						ZipEntry entry = zip.getEntry(entryName);
						Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
						CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
						
						for (AttenRelRef gmmRef : gmpeRefs) {
							String outputName = outputPrefix+"_"+gmmRef.name()+"_vs30_"+(int)vs30;
							if (inferred)
								outputName += "_inferred";
							else
								outputName += "_measured";
							outputName += ".csv";
							
							List<DiscretizedFunc> curves = new ArrayList<>();
							List<Double> weights = new ArrayList<>();
							double[] xVals = null;
							int colStart = -1;
							
							for (int i=0; i<tree.size(); i++) {
								LogicTreeBranch<?> branch = tree.getBranch(i);
								NGAW2_LogicTreeNode gmmNode = branch.requireValue(NGAW2_LogicTreeNode.class);
								if (gmmNode.getRef() == gmmRef) {
									weights.add(tree.getBranchWeight(i));
									
									if (xVals == null) {
										colStart = 3+branch.size();
										int numX = csv.getNumCols()-colStart;
										xVals = new double[numX];
										for (int j=0; j<xVals.length; j++)
											xVals[j] = csv.getDouble(0, colStart+j);
									}
									
									int row = i+1;
									double[] yVals = new double[xVals.length];
									for (int j=0; j<xVals.length; j++)
										yVals[j] = csv.getDouble(row, colStart+j);
									curves.add(new LightFixedXFunc(xVals, yVals));
								}
							}
							System.out.println("Loaded "+curves.size()+" curves for "+outputName);
							CSVFile<String> csvOut = TimJonEpistemicCurveCalc.buildCurvesCSV(percentiles, curves, weights);
							csvOut.writeToFile(new File(outputDir, outputName));
						}
					}
				}
			}
		}
		
		zip.close();
	}

}
