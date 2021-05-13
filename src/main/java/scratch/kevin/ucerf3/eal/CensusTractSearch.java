package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.param.Parameter;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;

public class CensusTractSearch {

	public static void main(String[] args) throws IOException {
		Location targetLoc = new Location(34.104587, -117.292812);
		
		Portfolio portfolio = Portfolio.createPortfolio(
				new File("/home/kevin/OpenSHA/UCERF3/eal/Porter-22-May-14-CA-ppty-90pct-Wills.txt"));
		
		int count = 0;
		
		HashSet<String> processed = new HashSet<String>();
		for (Asset asset : portfolio.getAssetList()) {
			String assetName = asset.getParameterList().getParameter(String.class, "AssetName").getValue();
			if (processed.contains(assetName))
				continue;
			processed.add(assetName);
			Location assetLoc = asset.getLocation();
			double dist = LocationUtils.horzDistanceFast(targetLoc, assetLoc);
			
			if (dist < 20d) {
				System.out.println("Asset "+assetName+", dist: "+dist);
//				for (Parameter<?> param : asset.getParameterList()) {
//					System.out.println("\t"+param.getName()+": "+param.getValue());
//				}
				count++;
			}
		}
		
		System.out.println("Found "+count+" matching");
	}

}
