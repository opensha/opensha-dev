package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.utils.FaultSectionDataWriter;

public class SubSectsWriterForRSQSim {

	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		NSHM23_DeformationModels[] dms = {
				NSHM23_DeformationModels.GEOLOGIC,
				NSHM23_DeformationModels.AVERAGE
		};
		File outputDir = new File("/tmp");
		
		String dateStr = new Date().toString();
		for (NSHM23_DeformationModels dm : dms) {
			List<? extends FaultSection> subSects = dm.build(fm);
			String prefix = fm.getFilePrefix()+"_"+dm.getFilePrefix()+"_sub_sections";
			GeoJSONFaultReader.writeFaultSections(new File(outputDir, prefix+".geojson"), subSects);
			List<String> metaData = new ArrayList<>();
			metaData.add("Fault Sub Sections file generated on "+dateStr+" by "+SubSectsWriterForRSQSim.class.getName());
			metaData.add("Fault Model: "+fm.getName());
			metaData.add("Deformation Model: "+dm.getName());
			FaultSectionDataWriter.writeSectionsToFile(subSects, metaData, new File(outputDir, prefix+".txt"), false);
			metaData.add("Note that upper seismogenic depths reflect aseismic reductions, and that slip rates " +
					"reported have a coupling coefficient applied and as such may be lower than the original " +
					"deformation model rate.");
			FaultSectionDataWriter.writeSectionsToFile(subSects, metaData, new File(outputDir, prefix+"_creep_reduced.txt"), true);
		}
	}

}
