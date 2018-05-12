package scratch.aftershockStatisticsETAS;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;

import com.google.common.io.Files;

public class GeoFeatureList extends ArrayList<GeoFeature> {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static boolean D = false; //debug
	
	public GeoFeatureList(){
		super();
	}

	public GeoFeatureList sortByMapLevel(){
		return this.sortByMapLevel(true);
	}
	public GeoFeatureList sortByMapLevel(boolean ascending){
		// sort by level
		this.sort(new Comparator<GeoFeature>(){
			@Override
			public int compare(GeoFeature a,GeoFeature b){
				if (ascending)
					return (a.mapLevel - b.mapLevel);
				else 
					return -(a.mapLevel - b.mapLevel); //sort in descending order
			}
		});
		return this;
	}

	public GeoFeatureList getFeaturesAboveLevel(int level){
		Iterator<GeoFeature> iter = this.iterator();
		while (iter.hasNext()){
			if (iter.next().mapLevel < level)
				iter.remove();
		}
		return this;
	}

	public GeoFeatureList getFeaturesInRegion(Region region){
		Iterator<GeoFeature> iter = this.iterator();
		while (iter.hasNext()){
			if (!iter.next().isInside(region))
				iter.remove();
		}
		return this;
	}

	/*
	 * Removes items from list that are within "distance" km of each other, preserving the item with 
	 * the higher mapLevel.
	 */
	public GeoFeatureList thinFeatures(double minDistance){
		double distance;
		int i = 0;
		while(i < this.size() ){
			int j = i + 1;
			while (j < this.size()){
				distance = this.get(i).distanceTo(this.get(j));
				if (distance < minDistance)
					this.remove(j);
				j++;
			}
			i++;	
		}
	return this;
	}

	public void addFeatureFromLine(String line){
		String[] elem = line.split(",");
		try{
			int mapLevel = Integer.parseInt(elem[1]);
			double lat = Double.parseDouble(elem[2]);
			double lon = Double.parseDouble(elem[3]);
			this.add(new GeoFeature("city", elem[0], new Location(lat, lon), mapLevel));
		
		} catch (Exception e){
			if(D){
				System.out.println("Unexpected line for GeoFeature string:");
				System.out.println(line);
				StringBuilder outline = new StringBuilder();
				for (String el : elem){
					outline.append(el + " ");
				}
				System.out.println(outline);
			}
		}
	}

	public static void main(String... args){
		GeoFeatureList cities = new GeoFeatureList();
		System.out.println("Loading Geographic information...");

		// load the data
		URL citiesURL = GeoFeatureList.class.getResource("worldcities1000.txt"); //updated CALIFORNIA added Mref column
		File cityFile = new File(citiesURL.getFile());
		
		List<String> lines = new ArrayList<String>();
		try{
			lines = Files.readLines(cityFile, Charset.defaultCharset());
		} catch (IOException e) {
			System.out.println("Couldn't load city information from " + cityFile);
		}

		//populate the feature list
		for (String line: lines){
//		for (int i = 0; i < 10; i++){
//			String line = lines.get(i);
			cities.addFeatureFromLine(line);
		}
		
		System.out.println(cities.size() + " cities added to list.");
		
		Region region = new Region(new Location(24,81), new Location(31,90));
		int minLevel = 1000;
		
		cities.getFeaturesInRegion(region);
		cities.getFeaturesAboveLevel(minLevel);
		cities.sortByMapLevel(false); //false --> descending order
		cities.thinFeatures(Math.sqrt(region.getExtent())/20);

		for(GeoFeature city : cities)
			System.out.println(city);
		
		System.out.println(cities.size() + " cities in mapped region.");
	}
}
