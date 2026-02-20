package scratch.kevin.pointSources.paperFigs2026;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.interp.DistanceInterpolator;
import org.opensha.sha.earthquake.PointSource.FocalMechRuptureSurfaceBuilder;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.DistanceDistributionCorrection;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.DistanceDistributionCorrection.FractileBin;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.DistanceDistributionCorrection.FractileDistances;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.DistanceDistributionCorrection.PrecomputedComparableDistances;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.util.FocalMech;

import com.google.common.base.Preconditions;

public class DistFileWriter {

	public static void main(String[] args) throws IOException {
		File dataDir = new File(ConstantsAndSettings.PAPER_DIR, "data");
		Preconditions.checkState(dataDir.exists() || dataDir.mkdir());
		
		boolean origProps = false;
		
		String propPrefix = origProps ? "original-rupture-properties" : "proposed-rupture-properties";
		
		// use these magnitudes
		EvenlyDiscretizedFunc magFunc = FaultSysTools.initEmptyMFD(3.01, 8.49); // this will still go to 8.55
		
		// use these distance bins
		DistanceInterpolator dists = DistanceInterpolator.get();
		
		FocalMechRuptureSurfaceBuilder surfBuilder = origProps ? ConstantsAndSettings.ORIG_SURF_BUILDER : ConstantsAndSettings.UPDATED_SURF_BUILDER;
		
		FocalMech[] mechs = {FocalMech.STRIKE_SLIP, FocalMech.REVERSE};
		
		DistanceDistributionCorrection corr = (DistanceDistributionCorrection)PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST.get();
		WeightedList<FractileBin> fractiles = corr.getFractiles();
		
		Location refLoc = new Location(0d, 0d);

		DecimalFormat oDF = new DecimalFormat("0.#");
		DecimalFormat magDF = new DecimalFormat("0.##");
		DecimalFormat lenDF = new DecimalFormat("0.0####");
		DecimalFormat dipDF = oDF;
		
		for (FocalMech mech : mechs) {
			List<String> header = new ArrayList<>();
			
			header.add("Magnitude");
			header.add("Length (km)");
			header.add("Down-Dip Width (km)");
			header.add("Dip (degrees)");
			header.add("Upper Depth (km)");
			header.add("Epicentral Distance (km)");
			
			boolean[] hws = mech == FocalMech.STRIKE_SLIP ? new boolean[] {false} : new boolean[] {true,false};
			
			if (mech != FocalMech.STRIKE_SLIP)
				header.add("Hanging Wall Fraction");
			for (boolean hw : hws) {
				String hwPrefix;
				if (mech == FocalMech.STRIKE_SLIP)
					hwPrefix = "";
				else if (hw)
					hwPrefix = "Hanging Wall ";
				else
					hwPrefix = "Footwall ";
				for (int f=0; f<fractiles.size(); f++) {
					FractileBin bin = fractiles.getValue(f);
					String binPrefix = hwPrefix+oDF.format(bin.minimum*100d)+"-"+oDF.format(bin.maximum*100d)+"%";
					header.add(binPrefix+" Rrup (km)");
					header.add(binPrefix+" Rjb (km)");
					header.add(binPrefix+" Rx (km)");
				}
			}
			
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine(header);
			
			for (int m=0; m<magFunc.size(); m++) {
				double mag = magFunc.getX(m);
				System.out.println("Doing "+mech+", M"+(float)mag);
				
				PointSurface surf = (PointSurface)surfBuilder.getSurface(refLoc, mag, mech, 0);
				
				List<String> commonHeader = List.of(
						magDF.format(mag),
						lenDF.format(surf.getAveLength()),
						lenDF.format(surf.getAveWidth()),
						dipDF.format(surf.getAveDip()),
						lenDF.format(surf.getAveRupTopDepth())
						);
				
				for (int d=0; d<dists.size(); d++) {
					double distance = dists.getDistance(d);
					
					FractileDistances fractileDists = corr.getFractileDistances(surf, distance);
					
					List<String> line = new ArrayList<>(header.size());
					line.addAll(commonHeader);
					line.add((float)distance+"");
					if (mech != FocalMech.STRIKE_SLIP)
						line.add((float)(1d-fractileDists.fractFootwall)+"");
					for (boolean hw : hws) {
						PrecomputedComparableDistances[] myDists = hw ? fractileDists.hangingWallDists : fractileDists.footwallDists;
						if (myDists == null) {
							Preconditions.checkState(!hw && fractileDists.fractFootwall == 0d,
									"Null dists for mech=%s, mag=%s, dist=%s, surf=%s, hw=%s, fractFW=%s",
									mech, mag, distance, surf, hw, fractileDists.fractFootwall);
							for (int i=0; i<fractiles.size(); i++) {
								line.add("NaN");
								line.add("NaN");
								line.add("NaN");
							}
							continue;
						}
						if (myDists.length == 1) {
							// all were the same
							PrecomputedComparableDistances[] copies = new PrecomputedComparableDistances[fractiles.size()]; 
							for (int i=0; i<copies.length; i++)
								copies[i] = myDists[0];
						} else {
							Preconditions.checkState(myDists.length == fractiles.size(),
									"Expected %s dists but have %s for mech=%s, mag=%s, dist=%s, surf=%s",
									fractiles.size(), myDists.length, mech, mag, distance, surf);
						}
						for (PrecomputedComparableDistances dist : myDists) {
							line.add((float)dist.getDistanceRup()+"");
							line.add((float)dist.getDistanceJB()+"");
							line.add((float)dist.getDistanceX()+"");
						}
					}
					csv.addLine(line);
				}
			}
			String outputPrefix = propPrefix + "-" + (mech == FocalMech.STRIKE_SLIP ? "strike-slip" : "dipping");
			csv.writeToFile(new File(dataDir, outputPrefix+".csv"));
		}
	}

}
