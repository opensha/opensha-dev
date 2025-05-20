package scratch.kevin.prvi25.figures;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import org.jfree.chart.ui.RectangleAnchor;
import org.opensha.commons.gui.plot.pdf.PDF_UTF8_FontMapper;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.RandomlySampledLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalRandomlySampledDeformationModelLevel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;

import com.google.common.base.Preconditions;
import com.itextpdf.awt.FontMapper;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.BaseColor;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfTemplate;
import com.itextpdf.text.pdf.PdfWriter;

import scratch.kevin.latex.LaTeXUtils;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class LogicTreeFigureWriter extends JPanel {

	public static void main(String[] args) throws IOException {
		Map<String, String> nameRemappings = new HashMap<>();
		nameRemappings.put(NSHM23_ScalingRelationships.WIDTH_LIMITED.getShortName(), "Wdth-Lmtd");
		nameRemappings.put(NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD.getShortName(), "W-L, CSD");
		nameRemappings.put(NSHM23_ScalingRelationships.LOGA_C4p2_SQRT_LEN.getShortName(), "Sqrt-Len");
		nameRemappings.put(PRVI25_CrustalRandomlySampledDeformationModelLevel.NAME, "Geologic Deformation Model Sample");
		
		File outputDir = new File(FIGURES_DIR, "logic_trees");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<LogicTree<?>> trees = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		
		LogicTree<?> crustalFaultTree = LogicTree.read(new File(CRUSTAL_DIR, "logic_tree.json"));
		trees.add(crustalFaultTree);
//		trees.set(trees.size()-1, reorderTree(trees.get(trees.size()-1), PRVI25_CrustalRandomlySampledDeformationModelLevel.NAME, 0));
		prefixes.add("crustal_inversion");
		LogicTree<?> crustalGriddedTree = new PRVI25_InvConfigFactory().getGridSourceTree(crustalFaultTree);
		trees.add(crustalGriddedTree);
		prefixes.add("crustal_gridded");
		
		// TODO switch back to reading it
//		LogicTree<LogicTreeNode> subFaultTree = LogicTree.read(new File(SUBDUCTION_DIR, "logic_tree.json"));
		LogicTree<LogicTreeNode> subFaultTree = LogicTree.buildExhaustive(PRVI25_LogicTreeBranch.levelsSubduction, true);
		trees.add(subFaultTree);
		prefixes.add("subduction_inversion");
		LogicTree<?> subGridTree = new PRVI25_InvConfigFactory().getGridSourceTree(subFaultTree);
		trees.add(subGridTree);
		prefixes.add("subduction_gridded");
		
		List<LogicTreeLevel<? extends LogicTreeNode>> crustalGMMLevels = PRVI25_LogicTreeBranch.levelsCrustalGMM;
		LogicTree<LogicTreeNode> crustalGMMTree = LogicTree.buildExhaustive(crustalGMMLevels, true);
		trees.add(crustalGMMTree);
		prefixes.add("gmm");
		
		List<LogicTreeLevel<? extends LogicTreeNode>> interfaceGMMLevels = PRVI25_LogicTreeBranch.levelsInterfaceGMM;
		LogicTree<LogicTreeNode> interfaceGMMTree = LogicTree.buildExhaustive(interfaceGMMLevels, true);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> slabGMMLevels = PRVI25_LogicTreeBranch.levelsSlabGMM;
		LogicTree<LogicTreeNode> slabGMMTree = LogicTree.buildExhaustive(slabGMMLevels, true);
		
		for (int i=0; i<trees.size(); i++) {
			LogicTree<?> tree = trees.get(i);
			String prefix = prefixes.get(i);
			System.out.println("Plotting "+prefix);
			
			LogicTreeFigureWriter panel = new LogicTreeFigureWriter(tree, false, nameRemappings);
			panel.write(outputDir, prefix, true, true);
			
			System.out.println();
		}
		
		FileWriter texFW = new FileWriter(new File(outputDir, "logic_tree_stats.tex"));
		
		texFW.write(LaTeXUtils.defineValueCommand("CrustalFaultBranches", LaTeXUtils.groupedIntNumber(crustalFaultTree.size()))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand("CrustalGriddedBranches", LaTeXUtils.groupedIntNumber(crustalGriddedTree.size()))+"\n");
		int crustalERFbranches = crustalFaultTree.size()*crustalGriddedTree.size();
		texFW.write(LaTeXUtils.defineValueCommand("CrustalCombinedERFBranches", LaTeXUtils.groupedIntNumber(crustalERFbranches))+"\n");
		int numWithoutRand = 1;
		HashSet<LogicTreeNode> crustalFaultNodesUsed = new HashSet<>();
		for (LogicTreeBranch<?> branch : crustalFaultTree)
			for (LogicTreeNode node : branch)
				crustalFaultNodesUsed.add(node);
		for (LogicTreeLevel<?> level : crustalFaultTree.getLevels()) {
			if (level instanceof RandomlySampledLevel<?>) {
				int numRand = level.getNodes().size();
				texFW.write(LaTeXUtils.defineValueCommand("CrustalFaultBranchesNumRand", LaTeXUtils.groupedIntNumber(numRand))+"\n");
			} else {
				int myNumNodes = 0;
				for (LogicTreeNode node : level.getNodes())
					if (crustalFaultNodesUsed.contains(node))
						myNumNodes++;
				Preconditions.checkState(myNumNodes > 0);
				numWithoutRand *= myNumNodes;
			}
		}
		texFW.write(LaTeXUtils.defineValueCommand("CrustalFaultBranchesWithoutRand", LaTeXUtils.groupedIntNumber(numWithoutRand))+"\n");
		int randMultiplier = crustalFaultTree.size() / numWithoutRand;
		texFW.write(LaTeXUtils.defineValueCommand("CrustalFaultBranchesRandMultiplier", LaTeXUtils.groupedIntNumber(randMultiplier))+"\n");
		
		texFW.write(LaTeXUtils.defineValueCommand("SubductionFaultBranches", LaTeXUtils.groupedIntNumber(subFaultTree.size()))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand("SubductionGriddedBranches", LaTeXUtils.groupedIntNumber(subGridTree.size()))+"\n");
		int subERFbranches = subFaultTree.size()*subGridTree.size();
		texFW.write(LaTeXUtils.defineValueCommand("SubductionCombinedERFBranches", LaTeXUtils.groupedIntNumber(subERFbranches))+"\n");
		
		int totalERFbranches = crustalERFbranches * subERFbranches;
		texFW.write(LaTeXUtils.defineValueCommand("CombinedERFBranches", LaTeXUtils.groupedIntNumber(totalERFbranches))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand("CombinedERFBranchesMillions", LaTeXUtils.groupedIntNumber(totalERFbranches*1e-6))+"\n");
		
		texFW.close();
	}
	
	private static LogicTree<LogicTreeNode> reorderTree(LogicTree<?> tree, String levelName, int targetIndex) {
		int foundIndex = -1;
		List<LogicTreeLevel<? extends LogicTreeNode>> newLevels = new ArrayList<>();
		for (int l=0; l<tree.getLevels().size(); l++) {
			LogicTreeLevel<?> level = tree.getLevels().get(l);
			if (level.getName().equals(levelName))
				foundIndex = l;
			newLevels.add(level);
		}
		Preconditions.checkState(foundIndex >= 0, "Didn't find level: %s", levelName);
		newLevels.add(targetIndex, newLevels.remove(foundIndex));
		List<LogicTreeBranch<LogicTreeNode>> newBranches = new ArrayList<>(tree.size());
		for (LogicTreeBranch<?> branch : tree) {
			LogicTreeBranch<LogicTreeNode> newBranch = new LogicTreeBranch<>(newLevels);
			for (LogicTreeNode node : branch)
				newBranch.setValue(node);
			newBranches.add(newBranch);
		}
		
		return LogicTree.fromExisting(newLevels, newBranches);
	}
	
	private static final int levelGap = 10;
	private static final int levelHGap = 10;
	private static final int choiceHGap = 5;
	private static final int lineGap = 5;
	private static final int minWidthPerNode = 100;
	private static final Font levelFont = new Font(Font.SANS_SERIF, Font.BOLD, 24);
	private static final Font choiceFont = new Font(Font.SANS_SERIF, Font.PLAIN, 18);
	private static final Font weightFont = new Font(Font.SANS_SERIF, Font.ITALIC, 18);
	private static final DecimalFormat weightDF = new DecimalFormat("0.###");
	private int levelFontHeight;
	private int choiceFontHeight;
	private int choiceWidth;
	private int lineHeight;
	private int levelHeight;
	private int width;
	private int height;
	private List<LogicTreeLevel<? extends LogicTreeNode>> includedLevels;
	private List<List<LogicTreeNode>> includedLevelChoices;

	private LogicTree<?> tree;
	private Map<String, String> nameRemappings;
	private double totalWeight;
	private Map<LogicTreeNode, Double> choiceWeights;
	private int maxNodes;
	
	public LogicTreeFigureWriter(LogicTree<?> tree, boolean includeSingleChoice, Map<String, String> nameRemappings) {
		this.tree = tree;
		if (nameRemappings == null)
			nameRemappings = Map.of();
		this.nameRemappings = nameRemappings;

		maxNodes = 1;
		int maxChoiceFontWidth = 0;
		int maxLevelFontWidth = 0;
		includedLevels = new ArrayList<>();
		includedLevelChoices = new ArrayList<>();
		totalWeight = 0d;
		choiceWeights = new HashMap<>();
		for (int i=0; i<tree.size(); i++)
			totalWeight += tree.getBranchWeight(i);
		for (LogicTreeLevel<?> level : tree.getLevels()) {
			// figure out how many choices I have
			List<LogicTreeNode> uniqueChoices = new ArrayList<>();
			HashSet<LogicTreeNode> uniqueChoicesSet = new HashSet<>();
			for (LogicTreeBranch<?> branch : tree) {
				LogicTreeNode node = branch.requireValue(level.getType());
				double weight = tree.getBranchWeight(branch);
				if (!uniqueChoicesSet.contains(node)) {
					uniqueChoices.add(node);
					uniqueChoicesSet.add(node);
					choiceWeights.put(node, weight);
				} else {
					choiceWeights.put(node, choiceWeights.get(node)+weight);
				}
			}
			System.out.println(level.getName()+" has "+uniqueChoices.size()+" unique nodes");
			if (includeSingleChoice || uniqueChoices.size() > 1) {
				includedLevels.add(level);
				includedLevelChoices.add(uniqueChoices);
				maxLevelFontWidth = Integer.max(maxLevelFontWidth, new JPanel().getFontMetrics(levelFont).stringWidth(remapped(level.getName())));
				if (!(level instanceof RandomlySampledLevel<?>)) {
					maxNodes = Integer.max(maxNodes, uniqueChoices.size());
					for (LogicTreeNode node : uniqueChoices) {
						String name = remapped(node.getShortName());
						FontMetrics metrics = new JPanel().getFontMetrics(choiceFont);
						maxChoiceFontWidth = Integer.max(maxChoiceFontWidth, metrics.stringWidth(name));
					}
				}
			}
		}
		
		System.out.println("We have "+includedLevelChoices.size()+" levels and at most "+maxNodes+" nodes");
		
		FontMetrics metrics = new JPanel().getFontMetrics(levelFont);
		levelFontHeight = metrics.getHeight();
		System.out.println("Level font height in pixels: " + levelFontHeight);
		metrics = new JPanel().getFontMetrics(choiceFont);
		choiceFontHeight = metrics.getHeight();
		System.out.println("Choice font height in pixels: " + choiceFontHeight);
		lineHeight = (int)(1.5d*levelFontHeight + 0.5);
		levelHeight = levelFontHeight + lineGap*2 + lineHeight + choiceFontHeight + choiceFontHeight;

		choiceWidth = Integer.max(minWidthPerNode, maxChoiceFontWidth+choiceHGap);
		width = Integer.max((int)(choiceWidth*maxNodes + 0.5*choiceWidth + 0.5),
				maxLevelFontWidth+2*levelHGap);
		height = levelHeight * includedLevelChoices.size() + levelGap * includedLevels.size();
		
		System.out.println("Calculated dimensions: "+width+"x"+height);
		
		setSize(width, height);
		setBackground(Color.WHITE);
	}
	
	private String remapped(String name) {
		if (nameRemappings.containsKey(name))
			return nameRemappings.get(name);
		return name;
	}
	
	@Override
	protected void paintComponent(Graphics g) {
		super.paintComponent(g);
		
		int y = (int)(0.5*levelGap + 0.5);
		int middleX = (int)(0.5*width + 0.5);
		
		Graphics2D g2d = (Graphics2D) g;
		float lineWidth = 3f;
		BasicStroke lineStroke = new BasicStroke(lineWidth);
		BasicStroke randLineStroke = new BasicStroke(lineWidth, BasicStroke.CAP_BUTT,
				BasicStroke.JOIN_BEVEL,0,new float[] {Float.min(6, Float.max(lineWidth*0.7f, 3))},0);
		g2d.setStroke(new BasicStroke(3));
		
		for (int l=0; l<includedLevels.size(); l++) {
			LogicTreeLevel<? extends LogicTreeNode> level = includedLevels.get(l);
			String name = remapped(level.getName());
			g2d.setFont(levelFont);
			drawText(g2d, name, middleX, y, RectangleAnchor.TOP);
			
			y += levelFontHeight;
			y += lineGap;
			List<LogicTreeNode> nodes = includedLevelChoices.get(l);
			int botLineY = y+lineHeight;
			int topChoiceY = botLineY + lineGap;
			int topWeightY = topChoiceY + choiceFontHeight;
			if (level instanceof RandomlySampledLevel<?> && nodes.size() > maxNodes && includedLevels.size() > 1) {
				// figure out how many branches we have without this
				HashSet<String> uniquesWithout = new HashSet<>();
				for (LogicTreeBranch<?> branch : tree) {
					String str = "";
					for (int l1=0; l1<branch.size(); l1++) {
						if (branch.getLevel(l1) == level)
							continue;
						str += "_"+branch.getValue(l1).getFilePrefix()+"_";
					}
					uniquesWithout.add(str);
//					System.out.println(str);
				}
				int prefNumLines;
				if (tree.size() % uniquesWithout.size() == 0) {
					int numPer = tree.size() / uniquesWithout.size();
					g2d.setFont(choiceFont);
					prefNumLines = numPer;
					drawText(g2d, numPer+" samples per branch, "+nodes.size()+" in total", middleX, topChoiceY, RectangleAnchor.TOP);
				} else {
					prefNumLines = nodes.size();
					drawText(g2d, nodes+" samples", middleX, topChoiceY, RectangleAnchor.TOP);
				}
				if (prefNumLines > maxNodes)
					prefNumLines = maxNodes;
				int myNodesWidth = choiceWidth*prefNumLines;
				int sideBuffer = (int)(0.5*(width - myNodesWidth));
				int leftX = sideBuffer;
				g2d.setStroke(randLineStroke);
				for (int i=0; i<prefNumLines; i++) {
					int rightX = leftX + choiceWidth;
					int choiceCenterX = (int)(0.5*(leftX + rightX));
					int topX = middleX + (int)(0.2*(choiceCenterX - middleX));
					g2d.drawLine(topX, y, choiceCenterX, botLineY);
					leftX = rightX;
				}
				g2d.setStroke(lineStroke);
				double minWeight = Double.POSITIVE_INFINITY;
				double maxWeight = Double.NEGATIVE_INFINITY;
				for (LogicTreeNode node : nodes) {
					double weight = node.getNodeWeight(null);
					minWeight = Math.min(minWeight, weight);
					maxWeight = Math.max(maxWeight, weight);
				}
				g2d.setFont(weightFont);
				if ((float)minWeight == (float)maxWeight)
					drawText(g2d, "(equally weighted)", middleX, topWeightY, RectangleAnchor.TOP);
				else
					drawText(g2d, "("+weightDF.format(minWeight)+"-"+weightDF.format(maxWeight)+")", middleX, topWeightY, RectangleAnchor.TOP);
			} else {
				int myNodesWidth = choiceWidth*nodes.size();
				int sideBuffer = (int)(0.5*(width - myNodesWidth));
				int leftX = sideBuffer;
				for (LogicTreeNode node : nodes) {
					int rightX = leftX + choiceWidth;
					int choiceCenterX = (int)(0.5*(leftX + rightX));
					int topX = middleX + (int)(0.2*(choiceCenterX - middleX));
					g2d.drawLine(topX, y, choiceCenterX, botLineY);
					leftX = rightX;
					g2d.setFont(choiceFont);
					drawText(g2d, remapped(node.getShortName()), choiceCenterX, topChoiceY, RectangleAnchor.TOP);
					double weight = choiceWeights.get(node)/totalWeight;
					g2d.setFont(weightFont);
					drawText(g2d, "("+weightDF.format(weight)+")", choiceCenterX, topWeightY, RectangleAnchor.TOP);
				}
			}
//			y += lineHeight;
//			y += lineGap;
//			y += choiceFontHeight;
//			y += choiceFontHeight;
			y = topWeightY + choiceFontHeight + levelGap;
		}
//		// Custom drawing code goes here
//		g2d.setStroke(new BasicStroke(2));
//		g2d.setColor(Color.BLACK);
//		g2d.drawLine(50, 50, 150, 150);
//		g2d.setFont(levelFont);
//		
//		drawText(g2d, "Test Text", 200, 0, RectangleAnchor.TOP);
		
	}
	
	private static void drawText(Graphics2D g2d, String text, int x, int y, RectangleAnchor anchor) {
		FontMetrics metrics = g2d.getFontMetrics();
		int textWidth = metrics.stringWidth(text);
		int textHeight = metrics.getHeight();
		int drawX = x;
		int drawY = y;

		switch (anchor) {
			case CENTER:
				drawX = x - textWidth / 2;
				drawY = y + metrics.getAscent() / 2 - textHeight / 4;
				break;
			case TOP_LEFT:
				drawX = x;
				drawY = y + metrics.getAscent();
				break;
			case TOP_RIGHT:
				drawX = x - textWidth;
				drawY = y + metrics.getAscent();
				break;
			case BOTTOM_LEFT:
				drawX = x;
				drawY = y - metrics.getDescent();
				break;
			case BOTTOM_RIGHT:
				drawX = x - textWidth;
				drawY = y - metrics.getDescent();
				break;
			case TOP:
				drawX = x - textWidth / 2;
				drawY = y + metrics.getAscent();
				break;
			case BOTTOM:
				drawX = x - textWidth / 2;
				drawY = y - metrics.getDescent();
				break;
			case LEFT:
				drawX = x;
				drawY = y + metrics.getAscent() / 2 - textHeight / 4;
				break;
			case RIGHT:
				drawX = x - textWidth;
				drawY = y + metrics.getAscent() / 2 - textHeight / 4;
				break;
		}

		g2d.drawString(text, drawX, drawY);
	}
	
	public void write(File outputDir, String prefix, boolean writePNG, boolean writePDF) throws IOException {
		if (writePNG) {
			BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			Graphics2D g2 = image.createGraphics();

			// Enable anti-aliasing for better quality
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

			// Render the panel onto the image
			this.setSize(width, height);
			this.paint(g2);
			g2.dispose();
			
			// Write the BufferedImage to a PNG file
			File outputFile = new File(outputDir, prefix+".png");
			ImageIO.write(image, "png", outputFile);
		}

		if (writePDF) {
			// step 1
			com.itextpdf.text.Rectangle pageSize = new com.itextpdf.text.Rectangle(width, height);
			pageSize.setBackgroundColor(BaseColor.WHITE);
			Document metadataDocument = new Document(pageSize);
			metadataDocument.addAuthor("OpenSHA");
			metadataDocument.addCreationDate();
			
			try {
				// step 2
				PdfWriter writer;

				writer = PdfWriter.getInstance(metadataDocument,
						new BufferedOutputStream(new FileOutputStream(new File(outputDir, prefix+".pdf"))));
				// step 3
				metadataDocument.open();
				// step 4
				PdfContentByte cb = writer.getDirectContent();
				
				PdfTemplate tp = cb.createTemplate(width, height);
				
				FontMapper fontMapper = new PDF_UTF8_FontMapper();
				PdfGraphics2D g2d = new PdfGraphics2D(tp, width, height, fontMapper);
				g2d.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_NORMALIZE);
				g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
				this.paint(g2d);
				g2d.dispose();
				cb.addTemplate(tp, 0, 0);
			}
			catch (DocumentException de) {
				throw ExceptionUtils.asRuntimeException(de);
			}
			// step 5
			metadataDocument.close();
		}
	}

}
