package scratch.kevin.faultArrays;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

public class FaultArrayCalc {
	
	public static int[] SAF_PARENTS = {
			97,		// Imperial
			170,	// Brawley (Seismic Zone) alt 1
			295,	// San Andreas (Coachella) rev
			284,	// San Andreas (San Gorgonio Pass-Garnet HIll)
			283,	// San Andreas (San Bernardino S)
			282,	// San Andreas (San Bernardino N)
			301,	// San Andreas (Mojave S)
			286,	// San Andreas (Mojave N)
			287,	// San Andreas (Big Bend)
	};
	
	public static int[] SJC_PARENTS = {
			28,		// San Jacinto (Superstition Mtn)
			99,		// San Jacinto (Borrego)
			101,	// San Jacinto (Coyote Creek)
			293,	// San Jacinto (Anza) rev
			401,	// San Jacinto (Stepovers Combined)
			289,	// San Jacinto (San Jacinto Valley)
			119,	// San Jacinto (San Bernardino)
	};
	
	public static int[] ELSINORE_PARENTS = {
			104,	// Laguna Salada
			103,	// Elsinore (Coyote Mountains)
			102,	// Elsinore (Julian)
			299,	// Elsinore (Temecula)
			402,	// Elsinore (Stepovers Combined)
			296,	// Elsinore (Glen Ivy) rev
			236,	// Whittier alt 1
	};
	
	private static final double trace_discr = 0.1;
	
	public static LocationList getSAF_LinearTrace(Map<Integer, FaultSectionPrefData> allParents) {
		return getLinearMultiParentTrace(allParents, SAF_PARENTS);
	}
	
	public static LocationList getSJC_LinearTrace(Map<Integer, FaultSectionPrefData> allParents) {
		return getLinearMultiParentTrace(allParents, SJC_PARENTS);
	}
	
	public static LocationList getElsinoreLinearTrace(Map<Integer, FaultSectionPrefData> allParents) {
		return getLinearMultiParentTrace(allParents, ELSINORE_PARENTS);
	}
	
	private static LocationList getLinearMultiParentTrace(Map<Integer, FaultSectionPrefData> allParents, int[] parentIDs) {
		LocationList ret = new LocationList();
		
		Location first = null;
		for (int parentID : parentIDs) {
			FaultSectionPrefData sect = allParents.get(parentID);
			FaultTrace sampledTrace = sect.getStirlingGriddedSurface(trace_discr, false, false).getEvenlyDiscritizedUpperEdge();
			double az = LocationUtils.azimuth(sampledTrace.first(), sampledTrace.last());
			if (az > 45 && az < 225)
				// reverse if it's NW to SE
				sampledTrace.reverse();
			for (Location loc : sampledTrace) {
				if (first == null) {
					first = loc;
				} else {
					double azFromFirst = LocationUtils.azimuth(first, loc);
					double azFromPrev = LocationUtils.azimuth(ret.last(), loc);
					double angleDiff = angleDiff(azFromFirst, azFromPrev);
					if (angleDiff > 90)
						// this is moving backwards, skip it
						continue;
				}
				ret.add(loc);
			}
		}
		
		return ret;
	}
	
	private static double angleDiff(double angle1, double angle2) {
		double angleDiff = Math.abs(angle1 - angle2);
		while (angleDiff > 270)
			angleDiff -= 360;
		return Math.abs(angleDiff);
	}
	
	public static List<Location> calcArrayLocations(LocationList linearTrace, double spacing) {
		double curDist = spacing-10;
		
		List<Location> arrayLocs = new ArrayList<>();
		
		for (int i=1; i<linearTrace.size(); i++) {
			Location prevLoc = linearTrace.get(i-1);
			Location loc = linearTrace.get(i);
			
			double thisDist = LocationUtils.horzDistanceFast(prevLoc, loc);
			if (thisDist + curDist >= spacing) {
				// it goes here
				Preconditions.checkState(thisDist < spacing);
				double distLeft = spacing - curDist;
				double azimuth = Math.toRadians(LocationUtils.azimuth(prevLoc, loc));
				arrayLocs.add(LocationUtils.location(prevLoc, azimuth, distLeft));
				curDist = thisDist - distLeft;
			} else {
				curDist += thisDist;
			}
		}
		
		return arrayLocs;
	}
	
	public static double calcFaultNormalAzimuth(LocationList linearTrace, Location loc, double bufferDistance) {
		List<Double> azimuths = new ArrayList<>();
		for (int i=0; i<linearTrace.size(); i++) {
			Location tLoc = linearTrace.get(i);
			if (LocationUtils.horzDistanceFast(loc, tLoc) <= bufferDistance) {
				double az = i > 0 ? LocationUtils.azimuth(linearTrace.get(i-1), tLoc) : LocationUtils.azimuth(tLoc, linearTrace.get(i+1));
				azimuths.add(az);
			}
		}
		return FaultUtils.getAngleAverage(azimuths)+90;
	}

}
