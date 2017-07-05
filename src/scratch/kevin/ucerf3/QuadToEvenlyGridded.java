package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.sha.earthquake.observedEarthquake.parsers.ngaWest.NGAWestParser;
import org.opensha.sha.faultSurface.ApproxEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GriddedSubsetSurface;
import org.opensha.sha.faultSurface.GriddedSurfaceImpl;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;

public class QuadToEvenlyGridded {

	public static EvenlyGriddedSurface getEvenlyGridded(EvenlyGriddedSurface quad, double approxSpacing) {
		return getEvenlyGridded(quad, approxSpacing, approxSpacing);
	}

	public static ArrayList<EvenlyGriddedSurface> getEvenlyGriddedPatches(
			EvenlyGriddedSurface quads, double approxSpacingX, double approxSpacingY) {
		Preconditions.checkNotNull(quads, "input surface cannot be null");
		Preconditions.checkArgument(quads.getNumCols() >= 2, "quads must have at least 2 columns");
		Preconditions.checkArgument(quads.getNumRows() == 2, "quads must have 2 rows");

		int windowCol = 0;

		// separate into individual quad patches
		ArrayList<GriddedSubsetSurface> quadPatches = new ArrayList<GriddedSubsetSurface>();
		while (windowCol <= (quads.getNumCols()-2)) {
			GriddedSubsetSurface quad = new GriddedSubsetSurface(2, 2, 0, windowCol, quads);
			quadPatches.add(quad);
			windowCol++;
		}

		// need to calculate nRows for all
		double[] ddws = new double[quadPatches.size()+1];
		int ddwCount = 0;
		for (GriddedSubsetSurface quad : quadPatches) {
			ddws[ddwCount++] = LocationUtils.linearDistance(quad.get(0, 0), quad.get(1, 0));
		}
		// now the rightmost one
		GriddedSubsetSurface lastQuad = quadPatches.get(quadPatches.size()-1);
		ddws[ddwCount++] = LocationUtils.linearDistance(lastQuad.get(0, 1), lastQuad.get(1, 1));
		double avgDDW = StatUtils.mean(ddws);
		int nRows = 1 + Math.round((float) (avgDDW / approxSpacingY));
		if (nRows < 2)
			nRows = 2;
		Preconditions.checkState(nRows >= 2);

		ArrayList<EvenlyGriddedSurface> gridPatches = new ArrayList<EvenlyGriddedSurface>();
		// go over each quad patch
		for (GriddedSubsetSurface quad : quadPatches) {
			Location topLeft =		quad.get(0, 0);
			Location bottomLeft =	quad.get(1, 0);
			Location topRight =		quad.get(0, 1);
			Location bottomRight =	quad.get(1, 1);

			// these are distances along the top and bottom edges of the quad, and the average
			double topLength = LocationUtils.linearDistance(topLeft, topRight);
			double bottomLength = LocationUtils.linearDistance(bottomLeft, bottomRight);
			double avgLength = 0.5 * (topLength + bottomLength);

			// calculate nCols using the average and approx spacing
			int nCols = 1 + Math.round((float) (avgLength / approxSpacingX));
			if (nCols < 2)
				nCols = 2;
//			Preconditions.checkState(nCols >= 2);

			// these are distances down the ends of the quads
			double leftAzumuth = LocationUtils.azimuth(topLeft, bottomLeft);
			double leftHorizDist = LocationUtils.horzDistance(topLeft, bottomLeft);
			double leftHorizDelta = leftHorizDist / (double)(nRows-1);
			double leftVertDist = LocationUtils.vertDistance(topLeft, bottomLeft);
			double leftVertDelta = leftVertDist / (double)(nRows-1);
			
			double rightAzumuth = LocationUtils.azimuth(topRight, bottomRight);
			double rightHorizDist = LocationUtils.horzDistance(topRight, bottomRight);
			double rightHorizDelta = rightHorizDist / (double)(nRows-1);
			double rightVertDist = LocationUtils.vertDistance(topRight, bottomRight);
			double rightVertDelta = rightVertDist / (double)(nRows-1);

			GriddedSurfaceImpl surf = new GriddedSurfaceImpl(nRows, nCols, approxSpacingX);
			surf.set(0, 0, topLeft);

			int lastRow = nRows-1;
			int lastCol = nCols-1;

			for (int row=0; row<nRows; row++) {
				// first fill in the first and last points in the row

				Location first, last;
				// set left and rightmost points
				if (row == 0) {
					first = topLeft;
					last = topRight;
				} else if (row == lastRow) {
					first = bottomLeft;
					last = bottomRight;
				} else {
					first = LocationUtils.location(surf.get(row-1, 0), new LocationVector(leftAzumuth, leftHorizDelta, leftVertDelta));
					last = LocationUtils.location(surf.get(row-1, lastCol), new LocationVector(rightAzumuth, rightHorizDelta, rightVertDelta));
				}

				surf.set(row, 0, first);
				surf.set(row, lastCol, last);

				double azimuth = LocationUtils.azimuth(first, last);
				double totVertDist = LocationUtils.vertDistance(first, last);
				double totHorizDist = LocationUtils.horzDistance(first, last);
				double verDistPerCol = totVertDist / (double)(nCols-1);
				double horizDistPerCol = totHorizDist / (double)(nCols-1);

				for (int col=1; col<(nCols-1); col++) {
					Location prevCol = surf.get(row, col-1);

					surf.set(row, col, LocationUtils.location(prevCol, new LocationVector(azimuth, horizDistPerCol, verDistPerCol)));
				}
			}

			gridPatches.add(surf);
		}

		return gridPatches;
	}

	public static EvenlyGriddedSurface stitchGriddedPatches(ArrayList<EvenlyGriddedSurface> gridPatches) {
		Preconditions.checkArgument(!gridPatches.isEmpty());
		// now patch them together
		if (gridPatches.size() == 1)
			return gridPatches.get(0);
				
		int nRows = -1;
		int totGriddedCols = 0;
		for (EvenlyGriddedSurface surf : gridPatches) {
			if (nRows < 0)
				nRows = surf.getNumRows();
			else
				Preconditions.checkArgument(nRows == surf.getNumRows());
			
			totGriddedCols += surf.getNumCols();
		}

		// remove duplicate points
		int adjustedNumCols = totGriddedCols - (gridPatches.size() - 1);

		EvenlyGriddedSurface surf = new GriddedSurfaceImpl(nRows, adjustedNumCols, gridPatches.get(0).getAveGridSpacing());

		int colOffset = 0;
		for (EvenlyGriddedSurface gridPatch : gridPatches) {
			for (int row=0; row<gridPatch.getNumRows(); row++)
				// leave the rightmost column out, will be set by the next patch
				for (int col=0; col<gridPatch.getNumCols()-1; col++)
					surf.set(row, col+colOffset, gridPatch.get(row, col));

			colOffset += gridPatch.getNumCols()-1;
		}

		// now fill in the last column
		EvenlyGriddedSurface lastPatch = gridPatches.get(gridPatches.size()-1);
		int lastCol = surf.getNumCols()-1;
		int lastPatchCol = lastPatch.getNumCols()-1;
		for (int row=0; row<nRows; row++) {
			surf.set(row, lastCol, lastPatch.get(row, lastPatchCol));
		}

		return surf;
	}

	public static EvenlyGriddedSurface getEvenlyGridded(EvenlyGriddedSurface quads, double approxSpacingX, double approxSpacingY) {
		ArrayList<EvenlyGriddedSurface> gridPatches = getEvenlyGriddedPatches(quads, approxSpacingX, approxSpacingY);
		
		return stitchGriddedPatches(gridPatches);
	}


	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File polTllDir = new File("src"+File.separator+"resources"+File.separator+"data"+File.separator+"ngaWest");
		int cnt = 0;
		ArrayList<EvenlyGriddedSurface> surfs = new ArrayList<EvenlyGriddedSurface>();
		for (RuptureSurface surf : NGAWestParser.loadPolTlls(polTllDir).values()) {
			if (surf instanceof EvenlyGriddedSurface) {
				EvenlyGriddedSurface gridded = (EvenlyGriddedSurface)surf;
				if (gridded.getNumCols() < 2 || gridded.getNumRows() < 2) {
					System.out.println("Weird size: "+gridded.getNumRows()+"x"+gridded.getNumCols());
					continue;
				}
				surfs.add(gridded);
			} else {
				for (RuptureSurface rupSurf : ((CompoundSurface)surf).getSurfaceList()) {
					EvenlyGriddedSurface gridded = (EvenlyGriddedSurface)rupSurf;
					if (gridded.getNumCols() < 2 || gridded.getNumRows() < 2) {
						System.out.println("Weird size: "+gridded.getNumRows()+"x"+gridded.getNumCols());
						continue;
					}
					surfs.add(gridded);
				}
			}
		}
		for (EvenlyGriddedSurface surf : surfs) {
			double az = LocationUtils.azimuth(surf.get(1, 1), surf.get(1, 0));
			surf.set(1, 0, LocationUtils.location(surf.get(1, 0), Math.toRadians(az), 20d));
//			EvenlyGriddedSurface gridded = surf;
			EvenlyGriddedSurface gridded = getEvenlyGridded(surf, 0.5d);
			
			
//			FaultTrace upper = new FaultTrace("upper");
//			FaultTrace lower = new FaultTrace("lower");
//			for (int i=0; i<surf.getNumCols(); i++) {
//				upper.add(surf.get(0, i));
//				lower.add(surf.get(surf.getNumRows()-1, i));
//			}
//			double az = LocationUtils.azimuth(lower.get(1), lower.get(0));
//			lower.set(0, LocationUtils.location(lower.get(0), Math.toRadians(az), 20d));
//			EvenlyGriddedSurface gridded = new ApproxEvenlyGriddedSurface(upper, lower, 0.5d);
			
			
			File file = new File("/tmp/surf"+(cnt++)+".txt");
			System.out.println("Writing "+file.getAbsolutePath());
			FileWriter fw = new FileWriter(file);
			for (int i=0; i<gridded.getNumRows(); i++) {
				for (int j=0; j<gridded.getNumCols(); j++) {
					Location loc = gridded.get(i, j);
					fw.write(loc.getLatitude()+"\t"+loc.getLongitude()+"\t"+loc.getDepth()+"\n");
				}
			}
			fw.close();
		}
	}

}
