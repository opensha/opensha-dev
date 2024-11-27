package scratch.kevin.nshm23;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class SectRIDistCSVWriter {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		BranchSectParticMFDs sectMFDs = sol.requireModule(BranchSectParticMFDs.class);
		
		double[] mags = {0d, 6d, 6.5d, 7d, 7.5d};
		
		double[] fractiles = {0.025,0.16, 0.84, 0.975};
		
		List<String> riHeader = List.of("Subsection Index", "Parent Section ID", "Subsection Name"
				,"Mean RI", "Min RI", "Max RI", "Std. Dev", "p2.5", "p16.0", "p84.0", "p97.5");
		List<String> rateHeader = new ArrayList<>(riHeader);
		for (int i=0; i<rateHeader.size(); i++)
			rateHeader.set(i, rateHeader.get(i).replace("RI", "Rate"));
		List<CSVFile<String>> riCSVs = new ArrayList<>(mags.length);
		List<CSVFile<String>> rateCSVs = new ArrayList<>(mags.length);
		for (int m=0; m<mags.length; m++) {
			CSVFile<String> riCSV = new CSVFile<>(true);
			riCSV.addLine(riHeader);
			riCSVs.add(riCSV);
			CSVFile<String> rateCSV = new CSVFile<>(true);
			rateCSV.addLine(rateHeader);
			rateCSVs.add(rateCSV);
		}
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.05, 8.45);
		List<String> mfdHeader = new ArrayList<>();
		mfdHeader.addAll(riHeader.subList(0, 3));
		for (int i=0; i<refMFD.size(); i++)
			mfdHeader.add((float)refMFD.getX(i)+"");
		CSVFile<String> mfdParticCSV = new CSVFile<>(true);
		mfdParticCSV.addLine(mfdHeader);
		CSVFile<String> mfdNuclCSV = new CSVFile<>(true);
		mfdNuclCSV.addLine(mfdHeader);
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			FaultSection sect = rupSet.getFaultSectionData(s);
			IncrementalMagFreqDist particMFD = sol.calcParticipationMFD_forSect(s, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			IncrementalMagFreqDist nuclMFD = sol.calcNucleationMFD_forSect(s, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			List<String> common = List.of(s+"", sect.getParentSectionId()+"", sect.getSectionName());
			List<String> partic = new ArrayList<>(common);
			List<String> nucl = new ArrayList<>(common);
			for (int i=0; i<refMFD.size(); i++) {
				partic.add((float)particMFD.getY(i)+"");
				nucl.add((float)nuclMFD.getY(i)+"");
			}
			mfdParticCSV.addLine(partic);
			mfdNuclCSV.addLine(nucl);
			
			for (boolean rate : new boolean[] {false,true}) {
				ArbDiscrEmpiricalDistFunc[] dists = new ArbDiscrEmpiricalDistFunc[mags.length];
				for (int m=0; m<dists.length; m++)
					dists[m] = new ArbDiscrEmpiricalDistFunc();
				
				for (int b=0; b<sectMFDs.getNumBranches(); b++) {
					double weight = sectMFDs.getBranchWeight(b);
					IncrementalMagFreqDist mfd = sectMFDs.getSectionMFD(b, s);
					double[] rates = new double[mags.length];
					for (Point2D pt : mfd)
						for (int m=0; m<mags.length; m++)
							if ((float)pt.getX() >= (float)mags[m])
								rates[m] += pt.getY();
					if (!rate)
						for (int m=0; m<mags.length; m++)
							rates[m] = 1d/rates[m];
					
					for (int m=0; m<mags.length; m++)
						dists[m].set(rates[m], weight);
				}
				for (int m=0; m<mags.length; m++) {
					List<String> line = new ArrayList<>(common);
					line.add((float)dists[m].getMean()+"");
					line.add((float)dists[m].getMinX()+"");
					line.add((float)dists[m].getMinY()+"");
					line.add((float)dists[m].getStdDev()+"");
					for (double fract : fractiles)
						line.add((float)dists[m].getInterpolatedFractile(fract)+"");
					if (rate)
						rateCSVs.get(m).addLine(line);
					else
						riCSVs.get(m).addLine(line);
				}
			}
		}
		
		File outputDir = new File("/tmp");
		String prefix = "nshm23_wus";
		mfdParticCSV.writeToFile(new File(outputDir, prefix+"_partic_mfds.csv"));
		mfdNuclCSV.writeToFile(new File(outputDir, prefix+"_nucl_mfds.csv"));
		for (int m=0; m<mags.length; m++) {
			String magPrefix;
			if (mags[m] == 0)
				magPrefix = prefix+"_supra_seis";
			else
				magPrefix = prefix+"_m"+(float)mags[m];
			rateCSVs.get(m).writeToFile(new File(outputDir, magPrefix+"_rate_dists.csv"));
			riCSVs.get(m).writeToFile(new File(outputDir, magPrefix+"_ri_dists.csv"));
		}
	}

}
