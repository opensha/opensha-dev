package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.magScalingRelations.MagAreaRelationship;
import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.MagScalingRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Ellsworth_B_WG02_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Shaw_2009_ModifiedMagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.CSVFile;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

public class ScalingRelCSVGen {

	public static void main(String[] args) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = new ArrayList<>();
		header.add("Fault ID");
		header.add("Fault Name");
		header.add("Length (km)");
		header.add("Upper Seismogenic Depth (km)");
		header.add("Lower Seismogenic Depth (km)");
		header.add("Dip");
		header.add("Width (km)");
		header.add("Subsection Length (km)");
		
		MagScalingRelationship[] scales = {
				new Ellsworth_B_WG02_MagAreaRel(),
				new HanksBakun2002_MagAreaRel(),
				new Shaw_2009_ModifiedMagAreaRel(),
				new WC1994_MagLengthRelationship()
		};
		for (MagScalingRelationship scale : scales) {
			header.add(scale.getName()+" 2-Subsection Min Mag");
			header.add(scale.getName()+" Full Segment Mag");
		}
		csv.addLine(header);
		
		for (FaultSection sect : NSHM23_FaultModels.WUS_FM_v1p4.getFaultSections()) {
			List<String> line = new ArrayList<>();
			
			line.add(sect.getSectionId()+"");
			line.add(sect.getSectionName());
			line.add((float)sect.getFaultTrace().getTraceLength()+"");
			line.add((float)sect.getOrigAveUpperDepth()+"");
			line.add((float)sect.getAveLowerDepth()+"");
			line.add((float)sect.getAveDip()+"");
			line.add((float)sect.getOrigDownDipWidth()+"");
			
			List<? extends FaultSection> subSects = sect.getSubSectionsList(0.5*sect.getOrigDownDipWidth(), 0, 2);
			double subSectLen = subSects.get(0).getFaultTrace().getTraceLength();
			line.add((float)subSectLen+"");
			
			double minMagArea = subSectLen*2d*sect.getOrigDownDipWidth();
			double maxMagArea = sect.getTraceLength()*sect.getOrigDownDipWidth();
			for (MagScalingRelationship scale : scales) {
				if (scale instanceof MagAreaRelationship) {
					line.add((float)scale.getMedianMag(minMagArea)+"");
					line.add((float)scale.getMedianMag(maxMagArea)+"");
				} else if (scale instanceof MagLengthRelationship) {
					line.add((float)scale.getMedianMag(subSectLen*2d)+"");
					line.add((float)scale.getMedianMag(sect.getTraceLength())+"");
				}
			}
			csv.addLine(line);
		}
		
		csv.writeToFile(new File("/tmp/nshm23_fm_v1p4_sect_scale_rels.csv"));
	}

}
