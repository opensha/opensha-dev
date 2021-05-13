package scratch.kevin.ucerf3;

import java.io.IOException;

import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;

public class TestBilinear {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(4, 4, -120, 34.0, 1.0);
		for (int xIndex=0; xIndex<4; xIndex++) {
			for (int yIndex=0; yIndex<4; yIndex++) {
				xyz.set(xIndex, yIndex, (double)(xIndex*yIndex));
			}
		}
		ArbDiscrGeoDataSet geo = new ArbDiscrGeoDataSet(true);
		
		for (double lat=34; lat<=37; lat+=0.1) {
			for (double lon=-120; lon<=-117; lon+=0.1) {
				geo.set(new Location(lat, lon), xyz.bilinearInterpolation(lon, lat));
			}
		}
		ArbDiscrGeoDataSet.writeXYZFile(geo, "/tmp/interp.txt");
	}

}
