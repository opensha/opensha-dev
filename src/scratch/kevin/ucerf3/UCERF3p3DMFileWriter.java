package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;

public class UCERF3p3DMFileWriter {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String date_str = "2013_04_09";
		String fm_str = "fm_3_1";
		File dmDir = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR.getParentFile(), "DeformationModels");
		File csvFile = new File(dmDir, "ucerf3p3_fm3p1_dms.csv");
		
		CSVFile<String> csv = CSVFile.readFile(csvFile, false);
		
		Map<Integer, DeformationSection> abm = Maps.newHashMap();
		Map<Integer, DeformationSection> neok = Maps.newHashMap();
		Map<Integer, DeformationSection> zeng = Maps.newHashMap();
		Map<Integer, DeformationSection> geol = Maps.newHashMap();
		Map<Integer, DeformationSection> geol_lower = Maps.newHashMap();
		Map<Integer, DeformationSection> geol_upper = Maps.newHashMap();
		
		Map<Integer, String> namesMap = Maps.newHashMap();
		
		// start at 2 because 2 header lines
		for (int i=2; i<csv.getNumRows(); i++) {
			List<String> line = csv.getLine(i);
			
			int cnt = 0;
			Integer parentID = Integer.parseInt(line.get(cnt++));
			String mini = line.get(cnt++);
			// correct format if necessary
			mini = DeformationModelFileParser.getMinisectionString(DeformationModelFileParser.parseMinisectionNumber(mini));
			
			String name = line.get(cnt++);
			namesMap.put(parentID, name);
			
			double lon1 = Double.parseDouble(line.get(cnt++));
			double lat1 = Double.parseDouble(line.get(cnt++));
			double lon2 = Double.parseDouble(line.get(cnt++));
			double lat2 = Double.parseDouble(line.get(cnt++));
			
			cnt++; // creep
			
			Location loc1 = new Location(lat1, lon1);
			Location loc2 = new Location(lat2, lon2);
			
			double abm_rate = Double.parseDouble(line.get(cnt++));
			double abm_rake = Double.parseDouble(line.get(cnt++));
			double neok_rate = Double.parseDouble(line.get(cnt++));
			double neok_rake = Double.parseDouble(line.get(cnt++));
			cnt++; // neok dip
			double zeng_rate = Double.parseDouble(line.get(cnt++));
			double zeng_rake = Double.parseDouble(line.get(cnt++));
			double geol_rate = Double.parseDouble(line.get(cnt++));
			double geol_lower_rate = Double.parseDouble(line.get(cnt++));
			double geol_upper_rate = Double.parseDouble(line.get(cnt++));
			double geol_rake = Double.parseDouble(line.get(cnt++));
			
			addLine(abm, parentID, loc1, loc2, abm_rate, abm_rake);
			addLine(neok, parentID, loc1, loc2, neok_rate, neok_rake);
			addLine(zeng, parentID, loc1, loc2, zeng_rate, zeng_rake);
			addLine(geol, parentID, loc1, loc2, geol_rate, geol_rake);
			addLine(geol_lower, parentID, loc1, loc2, geol_lower_rate, geol_rake);
			addLine(geol_upper, parentID, loc1, loc2, geol_upper_rate, geol_rake);
		}
		
		write(abm, new File(dmDir, "ABM_slip_rake_"+fm_str+"_"+date_str+".csv"), namesMap);
		write(neok, new File(dmDir, "neokinema_slip_rake_"+fm_str+"_"+date_str+".csv"), namesMap);
		write(zeng, new File(dmDir, "zeng_slip_rake_"+fm_str+"_b_bounded_"+date_str+".csv"), namesMap);
		write(geol, new File(dmDir, "geologic_slip_rake_"+fm_str+"_"+date_str+".csv"), namesMap);
		write(geol_lower, new File(dmDir, "geologic_slip_rake_"+fm_str+"_lowerbound_"+date_str+".csv"), namesMap);
		write(geol_upper, new File(dmDir, "geologic_slip_rake_"+fm_str+"_upperbound_"+date_str+".csv"), namesMap);
	}
	
	private static void addLine(Map<Integer, DeformationSection> dm, Integer parent,
			Location loc1, Location loc2, double slip, double rake) {
		DeformationSection sect = dm.get(parent);
		if (sect == null) {
			sect = new DeformationSection(parent);
			dm.put(parent, sect);
		}
		sect.add(loc1, loc2, slip, rake);
	}
	
	private static void write(Map<Integer, DeformationSection> dm, File file, final Map<Integer, String> namesMap) throws IOException {
		List<DeformationSection> sects = Lists.newArrayList(dm.values());
		Collections.sort(sects, new Comparator<DeformationSection>() {

			@Override
			public int compare(DeformationSection o1, DeformationSection o2) {
				return namesMap.get(o1.getId()).compareTo(namesMap.get(o2.getId()));
			}
		});
		
		DeformationModelFileParser.write(sects, file);
	}

}
