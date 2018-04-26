package scratch.aftershockStatisticsETAS;

import java.util.List;

import org.opensha.commons.util.cpt.CPT;

import wContour.Global.PolyLine;

public class ContourModel {
	private List<PolyLine> contours;
	private String name;
	private CPT cpt;

	public ContourModel(List<PolyLine> contours, String name, CPT cpt){
		this.contours = contours;
		this.name = name;
		this.cpt = cpt;
	}
	
	public List<PolyLine> getContours(){
		return contours;
	}
	
	public String getName(){
		return name;
	}
	
	public CPT getCPT(){
		return cpt;
	}
	
	public void setContours(List<PolyLine> contours){
		this.contours = contours;
	}

	public void setName(String name){
		this.name = name;
	}
	
	public void setCPT(CPT cpt){
		this.cpt = cpt;
	}
	
}


