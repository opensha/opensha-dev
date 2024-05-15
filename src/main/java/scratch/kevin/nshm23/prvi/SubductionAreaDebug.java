package scratch.kevin.nshm23.prvi;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultSection;

public class SubductionAreaDebug {

	public static void main(String[] args) throws IOException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/tmp/subduction_rup_set_FULL.zip"));
		
		List<? extends FaultSection> subSects = rupSet.getFaultSectionDataList();

		areaDebug(subSects.get(20));
		areaDebug(subSects.get(31));
		areaDebug(subSects.get(33));
		areaDebug(subSects.get(38));
		areaDebug(subSects.get(55));
		areaDebug(subSects.get(9));
	}
	
	private static void areaDebug(FaultSection sect) {
		System.out.println("Area debug for "+sect.getSectionId()+". "+sect.getSectionName());
		double sectArea = sect.getArea(false)*1e-6;
		System.out.println("\tFaultSection.getArea(false):\t"+(float)sectArea);
		for (double gridSpacing : new double[] {1d, 0.1d}) {
			EvenlyGriddedSurface surf = (EvenlyGriddedSurface) sect.getFaultSurface(gridSpacing, false, false);
			System.out.println("\tSurface area, gridSpacing:\t"+(float)gridSpacing);
			System.out.println("\t\tEvenlyGriddedSurface.getArea():\t"+surf.getArea()+" ("+pDiff(surf.getArea(), sectArea)+")");
			// now calculate our own
			double calcArea = 0d;
			int rows = surf.getNumRows();
			int cols = surf.getNumCols();
			for (int col=0; col<cols; col++) {
				Location top = surf.get(0, col);
				Location bot = surf.get(rows-1, col);
				double colLen = LocationUtils.linearDistanceFast(top, bot);
				double colWidth = 0d;
				if (col > 0) {
					double spacingToPrev = 0.5*(LocationUtils.linearDistanceFast(top, surf.get(0, col-1))
							+ LocationUtils.linearDistanceFast(bot, surf.get(rows-1, col-1)));
					colWidth += 0.5*spacingToPrev;
				}
				if (col < cols-1) {
					double spacingToNext = 0.5*(LocationUtils.linearDistanceFast(top, surf.get(0, col+1))
							+ LocationUtils.linearDistanceFast(bot, surf.get(rows-1, col+1)));
					colWidth += 0.5*spacingToNext;
				}
				calcArea += colLen*colWidth;
			}
			System.out.println("\t\tCalculated Area:\t"+(float)calcArea+" ("+pDiff(calcArea, sectArea)+")");
		}
	}
	
	private static final String pDiff(double test, double ref) {
		return pDF.format((test-ref)/ref);
	}
	
	private static DecimalFormat pDF = new DecimalFormat("0.00%");

}
