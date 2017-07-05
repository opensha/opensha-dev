package scratch.kevin.portfolioLEC;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;
import org.opensha.sra.gui.portfolioeal.PortfolioEALCalculatorController;
import org.opensha.sra.vulnerability.Vulnerability;

import com.google.common.base.Preconditions;

public class PortfolioDebug {

	public static void main(String[] args) throws IOException {
//		Portfolio portfolio = Portfolio.createPortfolio(new File("/home/kevin/OpenSHA/portfolio_lec/Porter-22-May-14-CA-CAS4-90pct-Wills.txt"));
		Portfolio portfolio = Portfolio.createPortfolio(new File("/tmp/Porter-17-Jun-14-CA-CAS4-90pct-Wills.csv"));
		System.out.println("Portolio has "+portfolio.getAssetList().size()+" assets");
		
		File vulnFile = new File("/home/kevin/OpenSHA/portfolio_lec/2014_05_16b_VUL06.txt");
		System.out.println("trying to load vulnerabilities from: "+vulnFile.getAbsolutePath());
		HashMap<String, Vulnerability> vulns = PortfolioEALCalculatorController.getVulnerabilities(vulnFile);
		System.out.println("DONE loading vulns.");
		
		for (Asset asset : portfolio.getAssetList())
			Preconditions.checkNotNull(vulns.get(asset.getVulnModelName()), "No vuln found for: "+asset.getVulnModelName());
	}

}
