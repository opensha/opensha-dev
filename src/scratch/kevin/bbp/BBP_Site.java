package scratch.kevin.bbp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.sha.earthquake.FocalMechanism;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class BBP_Site {
	
	private String name;
	private Location loc;
	private double vs30;
	private double loPassFreq;
	private double hiPassFreq;

	public BBP_Site(String name, Location loc, double vs30, double loPassFreq, double hiPassFreq) {
		this.name = name;
		this.loc = loc;
		this.vs30 = vs30;
		this.loPassFreq = loPassFreq;
		this.hiPassFreq = hiPassFreq;
	}

	public static void writeToFile(File file, Collection<BBP_Site> sites) throws IOException {
		FileWriter fw = new FileWriter(file);
		
		fw.write("#SLong    SLat     RSN   Vs30(m/s) LoPass_Freq(Hz) HiPass_Freq(Hz)\n");
		for (BBP_Site site : sites)
			fw.write(site.loc.getLongitude()+" "+site.loc.getLatitude()+" "+site.name+" "
					+site.vs30+" "+site.loPassFreq+" "+site.hiPassFreq+"\n");
		
		fw.close();
	}
	
	public static List<BBP_Site> readFile(File file) throws IOException {
		List<BBP_Site> sites = new ArrayList<>();
		
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("#"))
				continue;
			String[] split = line.split("\\s+");
			Preconditions.checkState(split.length == 6);
			Location loc = new Location(Double.parseDouble(split[1]), Double.parseDouble(split[0]));
			String name = split[2];
			double vs30 = Double.parseDouble(split[3]);
			double loPassFreq = Double.parseDouble(split[4]);
			double hiPassFreq = Double.parseDouble(split[5]);
			
			sites.add(new BBP_Site(name, loc, vs30, loPassFreq, hiPassFreq));
		}
		
		return sites;
	}

	public String getName() {
		return name;
	}

	public Location getLoc() {
		return loc;
	}

	public double getVs30() {
		return vs30;
	}

	public double getLoPassFreq() {
		return loPassFreq;
	}

	public double getHiPassFreq() {
		return hiPassFreq;
	}

}
