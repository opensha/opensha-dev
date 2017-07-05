package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.Lists;

import scratch.UCERF3.griddedSeismicity.FaultPolyMgr;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.utils.FaultSystemIO;

public class FaultPolyDataWriter {

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		InversionFaultSystemRupSet rupSet = FaultSystemIO.loadInvRupSet(new File("/home/kevin/OpenSHA/UCERF3/downsample_tests/"
				+ "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip"));
		FaultPolyMgr polyManager = rupSet.getInversionTargetMFDs().getGridSeisUtils().getPolyMgr();
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine("Sect Index", "Sect Name", "Parent Sect ID", "Num Poly Locations", "Latitude", "Longitude");
		for (FaultSectionPrefData sect : rupSet.getFaultSectionDataList()) {
			Region region = polyManager.getPoly(sect.getSectionId());
			LocationList border = region.getBorder();
			List<String> line = Lists.newArrayList(sect.getSectionId()+"", sect.getSectionName(),
					sect.getParentSectionId()+"", border.size()+"");
			for (Location loc : border) {
				line.add(loc.getLatitude()+"");
				line.add(loc.getLongitude()+"");
				csv.addLine(line);
				line = Lists.newArrayList("", "", "", "");
			}
		}
		csv.writeToFile(new File("/tmp/fm3_1_polys.csv"));
	}

}
