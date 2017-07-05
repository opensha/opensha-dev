package scratch.kevin.portfolioLEC;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;
import org.opensha.sra.gui.portfolioeal.PortfolioEALCalculatorController;
import org.opensha.sra.vulnerability.Vulnerability;

public class VulnPortfolioMismatchTest {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File portfolioFile = new File("/tmp/Porter-06-Feb-2012-CA99ptcPof.txt");
		File vulnFile = new File("/tmp/VUL06.txt");
		
		System.out.println("Parsing Portfolio: "+portfolioFile.getName());
		Portfolio port = Portfolio.createPortfolio(portfolioFile);
		System.out.println("DONE.");
		System.out.println("Parsing Vulnerabilities: "+vulnFile.getName());
		HashMap<String, Vulnerability> vulnMap = PortfolioEALCalculatorController.getVulnerabilities(vulnFile);
		System.out.println("DONE.");
		
		HashSet<String> alreadyDones = new HashSet<String>();
		
		
		System.out.println("Mismatches:");
		for (Asset asset : port.getAssetList()) {
			String vulnName = asset.getVulnModelName();
			if (alreadyDones.contains(vulnName))
				continue;
			alreadyDones.add(vulnName);
			if (!vulnMap.containsKey(vulnName))
				System.out.println("\t"+vulnName);
		}
	}

}
