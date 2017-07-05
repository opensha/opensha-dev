package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FaultUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.SimpleFaultData;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.utils.FaultSystemIO;

public class ScottTestRuptures {
	
	private static ArrayList<EvenlyGriddedSurface> getSurfaces(
			FaultSystemRupSet rupSet, int rupID, double gridSpacing) {
		List<Integer> sectIndexes = rupSet.getSectionsIndicesForRup(rupID);
		
		// bin them by common parents
		ArrayList<ArrayList<FaultSectionPrefData>> binnedByParent =
				new ArrayList<ArrayList<FaultSectionPrefData>>();
		binnedByParent.add(new ArrayList<FaultSectionPrefData>());
		
		for (int sectIndex : sectIndexes) {
			FaultSectionPrefData sect = rupSet.getFaultSectionData(sectIndex);
			ArrayList<FaultSectionPrefData> curList = binnedByParent.get(binnedByParent.size()-1);
			if (!curList.isEmpty() &&
					curList.get(curList.size()-1).getParentSectionId() != sect.getParentSectionId()) {
				// this means that we have a previous list, and it's parent is different.
				curList = new ArrayList<FaultSectionPrefData>();
				binnedByParent.add(curList);
			}
			curList.add(sect);
		}
		
		ArrayList<EvenlyGriddedSurface> surfs = new ArrayList<EvenlyGriddedSurface>();
		for (ArrayList<FaultSectionPrefData> sects : binnedByParent) {
			ArrayList<SimpleFaultData> sfds = new ArrayList<SimpleFaultData>();
			for (FaultSectionPrefData sect : sects)
				sfds.add(sect.getSimpleFaultData(false));
			SimpleFaultData sfd = SimpleFaultData.getCombinedSimpleFaultData(sfds);
			StirlingGriddedSurface surf = new StirlingGriddedSurface(sfds, gridSpacing);
			surfs.add(surf);
		}
		
		return surfs;
	}
	
	private static void writeRupture(FaultSystemRupSet rupSet, int rupID, File file, double gridSpacing)
			throws IOException {
//		ArrayList<EvenlyGriddedSurface> surfs = getSurfaces(rupSet, rupID, gridSpacing);
		List<Integer> sectIndexes = rupSet.getSectionsIndicesForRup(rupID);
		
		// bin them by common parents
		ArrayList<ArrayList<FaultSectionPrefData>> binnedByParent =
				new ArrayList<ArrayList<FaultSectionPrefData>>();
		binnedByParent.add(new ArrayList<FaultSectionPrefData>());
		
		for (int sectIndex : sectIndexes) {
			FaultSectionPrefData sect = rupSet.getFaultSectionData(sectIndex);
			ArrayList<FaultSectionPrefData> curList = binnedByParent.get(binnedByParent.size()-1);
			if (!curList.isEmpty() &&
					curList.get(curList.size()-1).getParentSectionId() != sect.getParentSectionId()) {
				// this means that we have a previous list, and it's parent is different.
				curList = new ArrayList<FaultSectionPrefData>();
				binnedByParent.add(curList);
			}
			curList.add(sect);
		}
		
		ArrayList<EvenlyGriddedSurface> surfs = new ArrayList<EvenlyGriddedSurface>();
		ArrayList<Double> rakes = new ArrayList<Double>();
		for (ArrayList<FaultSectionPrefData> sects : binnedByParent) {
			Preconditions.checkState(!sects.isEmpty());
			ArrayList<SimpleFaultData> sfds = new ArrayList<SimpleFaultData>();
			ArrayList<Double> parentRakes = new ArrayList<Double>();
			for (FaultSectionPrefData sect : sects) {
				sfds.add(sect.getSimpleFaultData(false));
				parentRakes.add(sect.getAveRake());
			}
			StirlingGriddedSurface surf = new StirlingGriddedSurface(sfds, gridSpacing);
			surfs.add(surf);
			double rake = FaultUtils.getAngleAverage(parentRakes);
			if (rake > 180)
				rake -= 360;
			rakes.add(rake);
		}
		
		FileWriter fw = new FileWriter(file);
		
		fw.write("Probability = 0\n");
		fw.write("Magnitude = "+rupSet.getMagForRup(rupID)+"\n");
		fw.write("GridSpacing = "+(float)gridSpacing+"\n");
		fw.write("Sections = "+surfs.size()+"\n");
		
		for (int i=0; i<surfs.size(); i++) {
			EvenlyGriddedSurface surf = surfs.get(i);
			fw.write("# Section "+i+"\n");
			fw.write("NumRows = "+surf.getNumRows()+"\n");
			fw.write("NumCols = "+surf.getNumCols()+"\n");
			fw.write("#   Lat         Lon         Depth      Rake    Dip     Strike\n");
			
			double rake = rakes.get(i);
			double dip = surf.getAveDip();
			double strike = surf.getAveStrike();
			
			for (Location loc : surf)
				fw.write(loc.getLatitude()+"\t"+loc.getLongitude()+"\t"+loc.getDepth()
						+"\t"+rake+"\t"+dip+"\t"+strike+"\n");
		}
		
		
	}

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 * @throws ZipException 
	 */
	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		// TODO Auto-generated method stub
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/preComputedData/FaultSystemRupSets/UCERF3_GEOLOGIC.zip"));
		double gridSpacing = 1d;
		writeRupture(rupSet, 137768, new File("/tmp/rupture_137768.txt"), gridSpacing);
		writeRupture(rupSet, 406047, new File("/tmp/rupture_406047.txt"), gridSpacing);
		writeRupture(rupSet, 292499, new File("/tmp/rupture_292499.txt"), gridSpacing);
		
	}

}
