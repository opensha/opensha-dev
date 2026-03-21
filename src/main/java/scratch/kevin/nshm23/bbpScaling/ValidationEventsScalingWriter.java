package scratch.kevin.nshm23.bbpScaling;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.data.CSVFile;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;

public class ValidationEventsScalingWriter {

	public static void main(String[] args) throws IOException {
		NSHM23_ScalingRelationships[] scales = {
				NSHM23_ScalingRelationships.AVERAGE,
				NSHM23_ScalingRelationships.LOGA_C4p1,
				NSHM23_ScalingRelationships.LOGA_C4p2,
				NSHM23_ScalingRelationships.LOGA_C4p3,
				NSHM23_ScalingRelationships.WIDTH_LIMITED,
				NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD,
				NSHM23_ScalingRelationships.LOGA_C4p2_SQRT_LEN
		};
		
		List<String> names = new ArrayList<>();
		List<Double> mags = new ArrayList<>();
		List<Double> lengths = new ArrayList<>();
		List<Double> widths = new ArrayList<>();
		
		names.add("ch_v14_2_2");
		mags.add(5.4);
		lengths.add(6.3);
		widths.add(4.6);
		
		names.add("m5_5_rv_socal");
		mags.add(5.5);
		lengths.add(5.62);
		widths.add(5.62);
		
		names.add("m6_2_ss_socal");
		mags.add(6.2);
		lengths.add(17.8);
		widths.add(8.9);
		
		names.add("m6_6_rv_socal");
		mags.add(6.6);
		lengths.add(28.2);
		widths.add(14.1);
		
		names.add("m6_6_ss_socal");
		mags.add(6.6);
		lengths.add(28.2);
		widths.add(14.1);
		
		names.add("nr_v14_02_1");
		mags.add(6.7);
		lengths.add(20d);
		widths.add(27d);
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = new ArrayList<>();
		header.add("Event");
		header.add("Original Magnitude");
		header.add("Original Length");
		header.add("Original DDW");
		for (NSHM23_ScalingRelationships scale : scales) {
			header.add(scale.getShortName()+" Magnitude");
			header.add(scale.getShortName()+" DDW for Original Magnitude");
		}
		csv.addLine(header);
		
		for (int i=0; i<names.size(); i++) {
			List<String> line = new ArrayList<>(header.size());
			double mag = mags.get(i);
			double length = lengths.get(i);
			double width = widths.get(i);
			line.add(names.get(i));
			line.add((float)mag+"");
			line.add((float)length+"");
			line.add((float)width+"");
			
			for (NSHM23_ScalingRelationships scale : scales) {
				double scaledMag = scale.getMag(length*width*1e6, length*1e3, width*1e3, width*1e3, Double.NaN);
				line.add((float)scaledMag+"");
				double scaledWidth = width;
				double testMag = scaledMag;
				int iter = 0;
				while (!Precision.equals(testMag, mag, 0.001)) {
					double deltaMag = mag - testMag;
					// width in km and mag aren't too far away in units
					scaledWidth += 0.1*deltaMag;
					testMag = scale.getMag(length*scaledWidth*1e6, length*1e3, scaledWidth*1e3, scaledWidth*1e3, Double.NaN);
//					System.out.println("iter "+iter+" testMag="+testMag+", scaledWidth="+scaledWidth+", origScaleMag="+scaledMag+", origMag="+mag+", origWidth="+width);
					iter++;
//					if (iter > 50)
//						System.exit(0);
				}
				line.add((float)scaledWidth+"");
			}
			csv.addLine(line);
		}
		csv.writeToFile(new File("/tmp/bbp_validation_events_nshm23_scaling.csv"));
	}

}
