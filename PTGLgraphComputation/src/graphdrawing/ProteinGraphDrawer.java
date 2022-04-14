/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package graphdrawing;

import io.FileParser;
import proteinstructure.SSE;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.commons.io.IOUtils;
import org.apache.xmlgraphics.java2d.GraphicContext;
import org.apache.xmlgraphics.java2d.ps.EPSDocumentGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import proteingraphs.FoldingGraph;
import io.IO;
import java.util.Map;
import resultcontainers.PTGLNotationFoldResult;
import proteingraphs.Position2D;
import proteingraphs.ProtGraph;
import proteingraphs.SSEGraph;
import settings.Settings;
import proteingraphs.SpatRel;
import proteinstructure.AminoAcid;
import proteinstructure.Chain;
import tools.DP;
import tools.PlccUtilities;

/**
 *
 * @author spirit
 */
public class ProteinGraphDrawer {
    
    private static final String AUTHOR_SSE_WARNING = "Created with author classification of SSEs. Results can't be compared to graphs of other proteins!";

    /**
     * Draws the SEQ legend for the graph at the given position. This legend is not suitable for other folding graph notations but SEQ, because their edges are different.
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position (x, y) where to start drawing
     * @param pl the page layout that defined the starting positions and borders of the elements on the canvas
     * @param g the graph to draw
     * @return the x coordinate in the image where the legend ends (which is the left margin + the legend width).
     * This can be used to determine the minimal width of the total image (it has to be at least this value).
     */
    public static Integer drawLegendSEQ(SVGGraphics2D ig2, Position2D startPos, PageLayout pl, SSEGraph g) {
        Boolean drawAll = Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend_always_all");
        ig2.setFont(pl.getLegendFont());
        FontMetrics fontMetrics = ig2.getFontMetrics();
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        Integer spacer = 10;
        Integer pixelPosX = startPos.x;
        Integer vertWidth = pl.getVertDiameter();
        Integer vertOffset = pl.getVertDiameter() / 4 * 3;
        String label;
        label = "[Edges: ";
        ig2.setPaint(Color.BLACK);
        ig2.drawString(label, pixelPosX, startPos.y);
        pixelPosX += fontMetrics.stringWidth(label) + spacer;
        label = "seq_distance";
        ig2.setPaint(Color.BLACK);
        ig2.drawString(label, pixelPosX, startPos.y);
        pixelPosX += fontMetrics.stringWidth(label) + spacer;
        label = "]";
        ig2.setPaint(Color.BLACK);
        ig2.drawString(label, pixelPosX, startPos.y);
        pixelPosX += fontMetrics.stringWidth(label) + spacer;
        if (g.numVertices() > 0) {
            label = " [Vertices: ";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeHelix() || drawAll) {
            ProteinGraphDrawer.drawSymbolAlphaHelix(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "helix";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeBetaStrand() || drawAll) {
            ProteinGraphDrawer.drawSymbolBetaStrand(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "strand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeLigand() || drawAll) {
            ProteinGraphDrawer.drawSymbolLigand(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "ligand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeOther() || drawAll) {
            ProteinGraphDrawer.drawSymbolOtherSSE(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "other";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.numVertices() > 0) {
            label = "]";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        return pixelPosX;
    }

    /**
     * Draws the ADJ folding graph image of this graph and returns the DrawResult.
     * If 'nonProteinGraph' is true, this graph is considered a custom (=non-protein) graph and the color coding for vertices and edges is NOT used.
     * In that case, the graph is drawn black and white and the labels for the N- and C-termini are NOT drawn.
     *
     * @param pnfr a folding graph notation result
     * @return the DrawResult. You can write this to a file or whatever.
     */
    private static DrawResult drawFoldingGraphADJG2D(PTGLNotationFoldResult pnfr) {
        FoldingGraph fg = pnfr.getFoldingGraph();
        SSEGraph pg = fg.getParent();
        Integer startVertexInParent = fg.getMinimalVertexIndexInParentGraph();
        Integer endVertexInParent = fg.getMaximalVertexIndexInParentGraph();
        Integer shiftBack = startVertexInParent;
        List<Integer> fgVertexPosInParent = fg.getVertexIndexListInParentGraph();
        List<Integer> parentVertexPosInFG = pg.computeParentGraphVertexPositionsInFoldingGraph(fgVertexPosInParent, pg.getSize());
        Integer numVerts = endVertexInParent - startVertexInParent + 1;
        Boolean bw = false;
        PageLayout pl = new PageLayout(numVerts);
        Position2D vertStart = pl.getVertStart();
        Integer lineHeight = pl.textLineHeight;
        SVGGraphics2D ig2;
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        ig2 = new SVGGraphics2D(document);
        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);
        Font font = pl.getStandardFont();
        Font fontBold = pl.getStandardFontBold();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();
        
        String proteinHeader = "ADJ " + pg.getGraphType() + " folding graph " + fg.getFoldingGraphFoldName() + " (# " + fg.getFoldingGraphNumber() + ") of PDB entry " + pg.getPdbid() + ", chain " + pg.getChainid() + "";
        proteinHeader += "\n" + Settings.getImagesCreatorVersionLine();
        if (FileParser.getWhichSseInfo().equals("author")){
            proteinHeader += "\n" + AUTHOR_SSE_WARNING;
        }
        
        String notation = "ADJ notation: '" + pnfr.adjNotation + "'";
        Integer stringHeight = fontMetrics.getAscent();
        String sseNumberSeq;
        String sseNumberFoldingGraph;
        String sseNumberProteinGraph;
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            DrawTools.drawStringLineBreaks(ig2, proteinHeader, pl.headerStart.x, pl.headerStart.y);
            if (Settings.getBoolean("PTGLgraphComputation_B_print_notations_on_fg_images")) {
                ig2.drawString(notation, pl.headerStart.x, pl.headerStart.y + lineHeight);
            }
        }
        Integer k;
        Integer l;
        Shape shape;
        Arc2D.Double arc;
        ig2.setStroke(new BasicStroke(2));
        Integer edgeType;
        Integer leftVert;
        Integer rightVert;
        Integer leftVertPosX;
        Integer rightVertPosX;
        Integer arcWidth;
        Integer arcHeight;
        Integer arcTopLeftX;
        Integer arcTopLeftY;
        Integer spacerX;
        Integer spacerY;
        for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
            for (Integer j = i + 1; j <= endVertexInParent; j++) {
                k = parentVertexPosInFG.get(i);
                l = parentVertexPosInFG.get(j);
                if (k < 0 || l < 0) {
                    continue;
                }
                if (fg.containsEdge(k, l)) {
                    edgeType = fg.getContactSpatRel(k, l);
                    if (edgeType.equals(SpatRel.PARALLEL)) {
                        ig2.setPaint(Color.RED);
                    } else if (edgeType.equals(SpatRel.ANTIPARALLEL)) {
                        ig2.setPaint(Color.BLUE);
                    } else if (edgeType.equals(SpatRel.MIXED)) {
                        ig2.setPaint(Color.GREEN);
                    } else if (edgeType.equals(SpatRel.LIGAND)) {
                        ig2.setPaint(Color.MAGENTA);
                    } else if (edgeType.equals(SpatRel.BACKBONE)) {
                        ig2.setPaint(Color.ORANGE);
                    } else {
                        ig2.setPaint(Color.LIGHT_GRAY);
                    }
                    if (bw) {
                        ig2.setPaint(Color.LIGHT_GRAY);
                    }
                    if (i < j) {
                        leftVert = i;
                        rightVert = j;
                    } else {
                        leftVert = j;
                        rightVert = i;
                    }
                    leftVertPosX = pl.getVertStart().x + ((leftVert - shiftBack) * pl.vertDist);
                    rightVertPosX = pl.getVertStart().x + ((rightVert - shiftBack) * pl.vertDist);
                    arcWidth = rightVertPosX - leftVertPosX;
                    arcHeight = arcWidth / 2;
                    arcTopLeftX = leftVertPosX;
                    arcTopLeftY = pl.getVertStart().y - arcHeight / 2;
                    spacerX = pl.vertRadius;
                    spacerY = 0;
                    arc = new Arc2D.Double(arcTopLeftX + spacerX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
                    shape = ig2.getStroke().createStrokedShape(arc);
                    ig2.fill(shape);
                }
            }
        }
        Ellipse2D.Double circle;
        Rectangle2D.Double rect;
        ig2.setStroke(new BasicStroke(2));
        for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
            if (pg.getVertex(i).isHelix()) {
                ig2.setPaint(Color.RED);
            } else if (pg.getVertex(i).isBetaStrand()) {
                ig2.setPaint(Color.BLACK);
            } else if (pg.getVertex(i).isLigandSSE()) {
                ig2.setPaint(Color.MAGENTA);
            } else if (pg.getVertex(i).isOtherSSE()) {
                ig2.setPaint(Color.GRAY);
            } else {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            Integer parentVertexPosInFoldingGraph = parentVertexPosInFG.get(i);
            if (parentVertexPosInFoldingGraph < 0) {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            if (bw) {
                ig2.setPaint(Color.GRAY);
            }
            SSE sse = pg.getVertex(i);
            if (sse.isBetaStrand()) {
                rect = new Rectangle2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(rect);
            } else if (sse.isLigandSSE()) {
                circle = new Ellipse2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.setStroke(new BasicStroke(3));
                ig2.draw(circle);
                ig2.setStroke(new BasicStroke(2));
            } else {
                circle = new Ellipse2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(circle);
            }
        }
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        if (!bw) {
            if (fg.getSize() > 0) {
                ig2.setFont(fontBold);
                ig2.drawString("N", vertStart.x - pl.vertDist, vertStart.y + 20);
                ig2.drawString("C", vertStart.x + numVerts * pl.vertDist, vertStart.y + 20);
                ig2.setFont(font);
            }
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {
            Integer printNth = 1;
            if (fg.getSize() > 99) {
                printNth = 2;
            }
            if (fg.getSize() > 999) {
                printNth = 3;
            }
            if (fg.getSize() > 0) {
                ig2.drawString("FG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                ig2.setPaint(Color.LIGHT_GRAY);
                ig2.drawString("SQ", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                ig2.drawString("PG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + lineHeight + (stringHeight / 4));
                ig2.setPaint(Color.BLACK);
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
                if ((i + 1) % printNth == 0) {
                    sseNumberFoldingGraph = "" + (parentVertexPosInFG.get(i) >= 0 ? (parentVertexPosInFG.get(i) + 1) : "");
                    sseNumberSeq = "" + (pg.getVertex(i).getSSESeqChainNum());
                    sseNumberProteinGraph = "" + (i + 1);
                    stringHeight = fontMetrics.getAscent();
                    ig2.setPaint(Color.BLACK);
                    ig2.drawString(sseNumberFoldingGraph, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (stringHeight / 4));
                    ig2.setPaint(Color.LIGHT_GRAY);
                    ig2.drawString(sseNumberSeq, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                    ig2.drawString(sseNumberProteinGraph, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));
                    ig2.setPaint(Color.BLACK);
                }
            }
            if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
                ProteinGraphDrawer.drawLegendDEFADJ(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4)), pl, fg, false);
            }
        }
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());
        DrawResult drawRes = new DrawResult(ig2, roi);
        return drawRes;
    }

    /**
     * Draws a protein graph in all formats, returns a list of written files.
     * @param baseFilePathNoExt the base file path where to put the image (without dot and file extension)
     * @param drawBlackAndWhite whether to omit colors, only useful for non-protein graphs
     * @param formats an array of type DrawTools.IMAGEFORMAT. Do not include SVG, this will always be drawn anyways
     * @param pg the graph to draw
     * @return a map of formats to the corresponding output files written to disk
     */
    public static HashMap<DrawTools.IMAGEFORMAT, String> drawProteinGraph(String baseFilePathNoExt, Boolean drawBlackAndWhite, DrawTools.IMAGEFORMAT[] formats, ProtGraph pg, Map<Integer, String> vertexMarkings, List<String> ignoreChains) {
        //Map<Integer, String> vertexMarkings = new HashMap<>();
        //for(Integer i : pg.getVertexIndexList()) {
        //    vertexMarkings.put(i, i + "");
        //}
        //System.out.println("####SETTING MARKING####");
        
        DrawResult drawRes = ProteinGraphDrawer.drawProteinGraphG2D(drawBlackAndWhite, pg, vertexMarkings, ignoreChains);
        String svgFilePath = baseFilePathNoExt + ".svg";
        HashMap<DrawTools.IMAGEFORMAT, String> resultFilesByFormat = new HashMap<DrawTools.IMAGEFORMAT, String>();
        try {
            DrawTools.writeG2dToSVGFile(svgFilePath, drawRes);
            resultFilesByFormat.put(DrawTools.IMAGEFORMAT.SVG, svgFilePath);
            resultFilesByFormat.putAll(DrawTools.convertSVGFileToOtherFormats(svgFilePath, baseFilePathNoExt, drawRes, formats));
        } catch (IOException ex) {
            DP.getInstance().e("Could not write protein graph file : '" + ex.getMessage() + "'.");
        }
        if (! (Settings.getBoolean("PTGLgraphComputation_B_silent") || Settings.getBoolean("PTGLgraphComputation_B_only_essential_output"))) {
            StringBuilder sb = new StringBuilder();
            sb.append("      Output protein graph files: ");
            for (DrawTools.IMAGEFORMAT format : resultFilesByFormat.keySet()) {
                sb.append("(").append(format.toString()).append(" => ").append(resultFilesByFormat.get(format)).append(") ");
            }
            System.out.println(sb.toString());
        }
        return resultFilesByFormat;
    }
    
    
    /**
     * Will draw anyything that implements IDrawableGraph, protein graph-style.
     * @param baseFilePathNoExt where to draw it, file without extension
     * @param formats the image formats to draw (extension gets appended based on this)
     * @param pg the graph to draw
     * @return a list of written files
     */
    public static HashMap<DrawTools.IMAGEFORMAT, String> drawDrawableGraph(String baseFilePathNoExt, DrawTools.IMAGEFORMAT[] formats, IDrawableGraph pg, Map<Integer, String> vertexMarkings) {
        DrawResult drawRes = ProteinGraphDrawer.drawDrawableGraphG2D(pg, vertexMarkings);
        String svgFilePath = baseFilePathNoExt + ".svg";
        HashMap<DrawTools.IMAGEFORMAT, String> resultFilesByFormat = new HashMap<DrawTools.IMAGEFORMAT, String>();
        try {
            DrawTools.writeG2dToSVGFile(svgFilePath, drawRes);
            resultFilesByFormat.put(DrawTools.IMAGEFORMAT.SVG, svgFilePath);
            resultFilesByFormat.putAll(DrawTools.convertSVGFileToOtherFormats(svgFilePath, baseFilePathNoExt, drawRes, formats));
        } catch (IOException ex) {
            DP.getInstance().e("Could not write drawable graph file : '" + ex.getMessage() + "'.");
        }
        if (!Settings.getBoolean("PTGLgraphComputation_B_silent")) {
            StringBuilder sb = new StringBuilder();
            sb.append("      Output drawable graph files: ");
            for (DrawTools.IMAGEFORMAT format : resultFilesByFormat.keySet()) {
                sb.append("(").append(format.toString()).append(" => ").append(resultFilesByFormat.get(format)).append(") ");
            }
            DP.getInstance().p(sb.toString());
        }
        return resultFilesByFormat;
    }

    /**
     * Draws a SEQ folding graph in all formats, returns a list of written files.
     * @param baseFilePathNoExt the base file path where to put the image (without dot and file extension)
     * @param drawBlackAndWhite whether to omit colors, only useful for non-protein graphs
     * @param formats an array of type DrawTools.IMAGEFORMAT. Do not include SVG, this will always be drawn anyways
     * @param pnfr
     * @return a map of formats to the corresponding output files written to disk
     */
    public static HashMap<DrawTools.IMAGEFORMAT, String> drawFoldingGraphSEQ(String baseFilePathNoExt, Boolean drawBlackAndWhite, DrawTools.IMAGEFORMAT[] formats, PTGLNotationFoldResult pnfr) {
        DrawResult drawRes = ProteinGraphDrawer.drawFoldingGraphSEQG2D(pnfr);
        String svgFilePath = baseFilePathNoExt + ".svg";
        HashMap<DrawTools.IMAGEFORMAT, String> resultFilesByFormat = new HashMap<DrawTools.IMAGEFORMAT, String>();
        try {
            DrawTools.writeG2dToSVGFile(svgFilePath, drawRes);
            resultFilesByFormat.put(DrawTools.IMAGEFORMAT.SVG, svgFilePath);
            resultFilesByFormat.putAll(DrawTools.convertSVGFileToOtherFormats(svgFilePath, baseFilePathNoExt, drawRes, formats));
        } catch (IOException ex) {
            DP.getInstance().e("Could not write SEQ folding graph file : '" + ex.getMessage() + "'.");
        }
        if (! (Settings.getBoolean("PTGLgraphComputation_B_silent") || Settings.getBoolean("PTGLgraphComputation_B_only_essential_output"))) {
            StringBuilder sb = new StringBuilder();
            sb.append("          Output SEQ folding graph files: ");
            for (DrawTools.IMAGEFORMAT format : resultFilesByFormat.keySet()) {
                sb.append("(").append(format.toString()).append(" => ").append(resultFilesByFormat.get(format)).append(") ");
            }
            System.out.println(sb.toString());
        }
        return resultFilesByFormat;
    }

    /**
     * Draws a RED folding graph in all formats, returns a list of written files.
     * @param baseFilePathNoExt the base file path where to put the image (without dot and file extension)
     * @param drawBlackAndWhite whether to omit colors, only useful for non-protein graphs
     * @param formats an array of type DrawTools.IMAGEFORMAT. Do not include SVG, this will always be drawn anyways
     * @param pnfr
     * @return a map of formats to the corresponding output files written to disk
     */
    public static HashMap<DrawTools.IMAGEFORMAT, String> drawFoldingGraphRED(String baseFilePathNoExt, Boolean drawBlackAndWhite, DrawTools.IMAGEFORMAT[] formats, PTGLNotationFoldResult pnfr) {
        DrawResult drawRes = ProteinGraphDrawer.drawFoldingGraphREDG2D(pnfr);
        String svgFilePath = baseFilePathNoExt + ".svg";
        HashMap<DrawTools.IMAGEFORMAT, String> resultFilesByFormat = new HashMap<DrawTools.IMAGEFORMAT, String>();
        try {
            DrawTools.writeG2dToSVGFile(svgFilePath, drawRes);
            resultFilesByFormat.put(DrawTools.IMAGEFORMAT.SVG, svgFilePath);
            resultFilesByFormat.putAll(DrawTools.convertSVGFileToOtherFormats(svgFilePath, baseFilePathNoExt, drawRes, formats));
        } catch (IOException ex) {
            DP.getInstance().e("Could not write RED folding graph file : '" + ex.getMessage() + "'.");
        }
        if (! (Settings.getBoolean("PTGLgraphComputation_B_silent") || Settings.getBoolean("PTGLgraphComputation_B_only_essential_output"))) {
            StringBuilder sb = new StringBuilder();
            sb.append("          Output RED folding graph files: ");
            for (DrawTools.IMAGEFORMAT format : resultFilesByFormat.keySet()) {
                sb.append("(").append(format.toString()).append(" => ").append(resultFilesByFormat.get(format)).append(") ");
            }
            System.out.println(sb.toString());
        }
        return resultFilesByFormat;
    }

    /**
     * Draws the KEY legend for the graph at the given position. This legend is not suitable for other folding graph notations but KEY, because their edges and SSE symbols are different.
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position (x, y) where to start drawing
     * @return the x coordinate in the image where the legend ends (which is the left margin + the legend width).
     * This can be used to determine the minimal width of the total image (it has to be at least this value).
     */
    public static Integer drawLegendKEY(SVGGraphics2D ig2, Position2D startPos, PageLayout pl, SSEGraph g) {
        Boolean drawAll = Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend_always_all");
        ig2.setFont(pl.getLegendFont());
        FontMetrics fontMetrics = ig2.getFontMetrics();
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        Integer spacer = 10;
        Integer pixelPosX = startPos.x;
        Integer vertWidth = pl.getVertDiameter();
        Integer vertOffset = pl.getVertDiameter() / 4 * 3;
        String label;
        label = "[Edges: ";
        ig2.setPaint(Color.BLACK);
        ig2.drawString(label, pixelPosX, startPos.y);
        pixelPosX += fontMetrics.stringWidth(label) + spacer;
        label = "sequential";
        ig2.setPaint(Color.BLACK);
        ig2.drawString(label, pixelPosX, startPos.y);
        pixelPosX += fontMetrics.stringWidth(label) + spacer;
        label = "]";
        ig2.setPaint(Color.BLACK);
        ig2.drawString(label, pixelPosX, startPos.y);
        pixelPosX += fontMetrics.stringWidth(label) + spacer;
        ig2.setPaint(Color.BLACK);
        if (g.numVertices() > 0) {
            label = " [Vertices: ";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeHelix() || drawAll) {
            label = "helix";
            ig2.setPaint(Color.RED);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeBetaStrand() || drawAll) {
            label = "strand";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeLigand() || drawAll) {
            ig2.setPaint(Color.MAGENTA);
            label = "ligand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeOther() || drawAll) {
            ig2.setPaint(Color.GRAY);
            label = "other";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.numVertices() > 0) {
            label = "]";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        return pixelPosX;
    }

    /**
     * Draws the KEY folding graph image of this graph and returns the DrawResult.
     * If 'nonProteinGraph' is true, this graph is considered a custom (=non-protein) graph and the color coding for vertices and edges is NOT used.
     * In that case, the graph is drawn black and white and the labels for the N- and C-termini are NOT drawn.
     *
     * See http://ptgl.uni-frankfurt.de/cgi-bin/showpict.pl?topology=a&rep=3&protlist=7timA+7timB+&nmrlist=
     * and http://ptgl.uni-frankfurt.de/ptglhelp.html#key
     *
     * See the KEY beta FG #3 of 1GOS chain A as an example.
     *
     * @param pnfr a folding graph notation result
     * @return the DrawResult. You can write this to a file or whatever.
     */
    @Deprecated
    private static DrawResult drawFoldingGraphKEYG2DOld(PTGLNotationFoldResult pnfr) {
        FoldingGraph fg = pnfr.getFoldingGraph();
        SSEGraph pg = fg.getParent();
        Boolean debugVerbose = false;
        if (fg.isBifurcated()) {
            DP.getInstance().e("SSEGraph", "drawFoldingGraphKEYG2D: This FG is bifurcated, KEY notation not supported.");
            return null;
        }
        Integer leftMostVertexInParent = fg.getMinimalVertexIndexInParentGraph();
        Integer rightMostVertexInParent = fg.getMaximalVertexIndexInParentGraph();
        List<Integer> keyposParentIndicesSeqOrder = pnfr.keypos;
        Integer keystartFGIndex = pnfr.keyStartFG;
        if (keystartFGIndex == null) {
            DP.getInstance().e("SSEGraph", "drawFoldingGraphKEYG2D: keystartFGIndex is null");
        }
        fg.computeSpatialVertexOrdering();
        List<Integer> keyposFGIndicesSpatOrder = fg.getSpatOrder();
        if (keyposFGIndicesSpatOrder.size() != fg.getSize()) {
            if (fg.isASingleCycle()) {
                keyposFGIndicesSpatOrder = fg.getSpatialOrderingOfVertexIndicesForSingleCycleFG(keystartFGIndex);
                fg.setSpatOrder(keyposFGIndicesSpatOrder);
                if (keyposFGIndicesSpatOrder.size() != fg.getSize()) {
                    DP.getInstance().e("SSEGraph", "Could not draw KEY notation: Spatorder size " + keyposFGIndicesSpatOrder.size() + " of chain " + fg.getChainid() + " gt " + fg.getGraphType() + " FG " + fg.getFoldingGraphFoldName() + " (#" + fg.getFoldingGraphNumber() + ") does not match FG size " + fg.getSize() + " even when allowing circles. Parent verts of FG: " + IO.intListToString(fg.getVertexIndexListInParentGraph()) + ". KEY='" + pnfr.keyNotation + "'. " + (fg.isASingleCycle() ? "FG is a single cycle." : "FG is NOT a single cycle."));
                    return null;
                } else {
                    //System.err.println("Single cycle handling fixed it.");
                }
            } else {
                DP.getInstance().e("SSEGraph", "Could not draw KEY notation: Spatorder size " + keyposFGIndicesSpatOrder.size() + " of chain " + fg.getChainid() + " gt " + fg.getGraphType() + " FG " + fg.getFoldingGraphFoldName() + " (#" + fg.getFoldingGraphNumber() + ") does not match FG size " + fg.getSize() + ". Parent verts of FG: " + IO.intListToString(fg.getVertexIndexListInParentGraph()) + ". KEY='" + pnfr.keyNotation + "'. " + (fg.isASingleCycle() ? "FG is a single cycle." : "FG is NOT a single cycle."));
                return null;
            }
        }
        if (debugVerbose) {
            DP.getInstance().d("SSEGraph", fg.getQuickIDString() + " Spatorder size " + keyposFGIndicesSpatOrder.size() + " of chain " + fg.getChainid() + " gt " + fg.getGraphType() + " FG " + fg.getFoldingGraphFoldName() + " (#" + fg.getFoldingGraphNumber() + "). FG size " + fg.getSize() + ". Parent verts of FG: " + IO.intListToString(fg.getVertexIndexListInParentGraph()) + ". KEY='" + pnfr.keyNotation + "'.");
        }
        List<Integer> spatOrderseqDistToPrevious = new ArrayList<>();
        spatOrderseqDistToPrevious.add(null);
        Integer cur;
        Integer last;
        Integer spatPosCur;
        Integer spatPosLast;
        for (Integer i = 1; i < fg.getSize(); i++) {
            cur = i;
            last = i - 1;
            spatPosCur = keyposFGIndicesSpatOrder.indexOf(cur);
            spatPosLast = keyposFGIndicesSpatOrder.indexOf(last);
            spatOrderseqDistToPrevious.add(spatPosCur - spatPosLast);
        }
        List<Integer> fgVertexPosInParent = fg.getVertexIndexListInParentGraph();
        List<Integer> parentVertexPosInFG = ProtGraph.computeParentGraphVertexPositionsInFoldingGraph(fgVertexPosInParent, pg.getSize());
        Integer numVerts = fg.getSize();
        Boolean bw = false;
        PageLayout pl = new PageLayout(numVerts);
        pl.isForKEY = true;
        Position2D vertStart = pl.getVertStart();
        Integer lineHeight = pl.textLineHeight;
        Integer vertDist = 50;
        Integer vertHeight = 80;
        Integer vertWidth = 40;
        Integer vertStartX = pl.getVertStart().x;
        Integer vertStartY = pl.getVertStart().y;
        SVGGraphics2D ig2;
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        ig2 = new SVGGraphics2D(document);
        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);
        Font font = pl.getStandardFont();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();
        String proteinHeader = "KEY " + pg.getGraphType() + " folding graph " + fg.getFoldingGraphFoldName() + " (# " + fg.getFoldingGraphNumber() + ") of PDB entry " + pg.getPdbid() + ", chain " + pg.getChainid() + "";
        String notation = "KEY notation: '" + pnfr.keyNotation + "'";
        Integer stringHeight = fontMetrics.getAscent();
        String sseNumberSeq;
        String sseNumberFoldingGraph;
        String sseNumberProteinGraph;
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            ig2.drawString(proteinHeader, pl.headerStart.x, pl.headerStart.y);
            ig2.drawString(notation, pl.headerStart.x, pl.headerStart.y + lineHeight);
        }
        Integer[] shiftBack = new Integer[pg.getSize()];
        Integer shift = 0;
        for (int i = 0; i < parentVertexPosInFG.size(); i++) {
            if (parentVertexPosInFG.get(i) < 0) {
                shift++;
            }
            shiftBack[i] = shift;
        }
        boolean debugDrawingKEY = false;
        //if(fg.toShortString().equals("1gos-B-alpha-FG0[4V,3E]")) {
        //    debugDrawingKEY = true;
        //}
        Integer[] orientationsSpatOrder = new Integer[fg.getSize()];
        Integer[] orientationsSeqOrder = new Integer[fg.getSize()];
        List<Integer> keyposFGIndices = new ArrayList<>();
        Integer v;
        for (int i = 0; i < keyposParentIndicesSeqOrder.size(); i++) {
            v = parentVertexPosInFG.get(keyposParentIndicesSeqOrder.get(i));
            if (v < 0) {
                DP.getInstance().e("SSEGRaph", "Keypos parent index #" + i + " is " + v + ", skipping.");
            } else {
                keyposFGIndices.add(v);
            }
        }
        if (keyposFGIndices.size() != fg.getSize()) {
            DP.getInstance().e("SSEGRaph", "keyposFGIndices.size()=" + keyposFGIndices.size() + ", fg.size=" + fg.getSize() + ".");
        }
        Integer firstVertexSpatFGIndex = fg.getSpatOrder().get(0);
        orientationsSpatOrder[0] = FoldingGraph.ORIENTATION_UPWARDS;
        orientationsSeqOrder[firstVertexSpatFGIndex] = FoldingGraph.ORIENTATION_UPWARDS;
        if (debugDrawingKEY) {
            System.out.println("Setting orientationsSeqOrder[" + firstVertexSpatFGIndex + "] to " + FoldingGraph.ORIENTATION_UPWARDS + ".");
        }
        Integer currentVert;
        Integer lastVert;
        StringBuilder KEYNotation = new StringBuilder();
        String bracketStart = "{";
        String bracketEnd = "}";
        KEYNotation.append(bracketStart);
        KEYNotation.append(fg.getVertex(keystartFGIndex).getLinearNotationLabel());
        if (keyposFGIndices.size() > 1) {
            for (int i = 1; i < keyposFGIndicesSpatOrder.size(); i++) {
                if (i < (keyposFGIndicesSpatOrder.size())) {
                    KEYNotation.append(",");
                }
                currentVert = keyposFGIndicesSpatOrder.get(i);
                lastVert = keyposFGIndicesSpatOrder.get(i - 1);
                KEYNotation.append(currentVert - lastVert);
                Integer spatRel = fg.getContactSpatRel(lastVert, currentVert);
                if (Objects.equals(spatRel, SpatRel.PARALLEL)) {
                    KEYNotation.append("x");
                    orientationsSpatOrder[i] = orientationsSpatOrder[i - 1];
                    orientationsSeqOrder[currentVert] = orientationsSpatOrder[i];
                } else {
                    orientationsSpatOrder[i] = (Objects.equals(orientationsSpatOrder[i - 1], FoldingGraph.ORIENTATION_UPWARDS) ? FoldingGraph.ORIENTATION_DOWNWARDS : FoldingGraph.ORIENTATION_UPWARDS);
                    orientationsSeqOrder[currentVert] = orientationsSpatOrder[i];
                }
                KEYNotation.append(fg.getVertex(currentVert).getLinearNotationLabel());
            }
        }
        KEYNotation.append(bracketEnd);
        if (debugDrawingKEY) {
            DP.getInstance().d("orientationsSpatOrder=" + IO.integerArrayToString(orientationsSpatOrder));
            DP.getInstance().d("orientationsSeqOrder=" + IO.integerArrayToString(orientationsSeqOrder));
        }
        Integer edgeType;
        Integer connCenterX;
        Integer connCenterY;
        Integer leftVertSeq;
        Integer rightVertSeq;
        Integer leftVertSpat;
        Integer rightVertSpat;
        Integer leftVertPosX;
        Integer leftVertPosY;
        Integer rightVertPosX;
        Integer rightVertPosY;
        Integer connWidth;
        Integer connHeight;
        Integer connTopLeftX;
        Integer connTopLeftY;
        Integer spacerX;
        Integer spacerY;
        Integer iSpatIndex;
        Integer jSpatIndex;
        Boolean startUpwards;
        Integer relDrawDistToLast;
        ig2.setStroke(new BasicStroke(1));
        if (debugDrawingKEY) {
            System.out.println("orientationsSeqOrder=" + IO.integerArrayToString(orientationsSeqOrder) + ", orientationsSpatOrder=" + IO.integerArrayToString(orientationsSpatOrder));
        }
        int numContactsDrawn = 0;
        if (keyposFGIndices.size() > 1) {
            for (int i = 1; i < keyposFGIndicesSpatOrder.size(); i++) {
                currentVert = keyposFGIndicesSpatOrder.get(i);
                lastVert = keyposFGIndicesSpatOrder.get(i - 1);
                Integer spatRel = fg.getContactSpatRel(lastVert, currentVert);
                relDrawDistToLast = currentVert - lastVert;
                if (debugDrawingKEY) {
                    System.out.println("At i=" + 1 + ", currentVert=" + currentVert + ", lastVert=" + lastVert + "");
                }
                if (spatRel.equals(SpatRel.NONE)) {
                    if (debugDrawingKEY) {
                        System.out.println("  skipping");
                    }
                    continue;
                }
                if (debugDrawingKEY) {
                    System.out.println("  spatRel=" + SpatRel.getString(spatRel) + ", drawing edge " + lastVert + " => " + currentVert + "...");
                }
                boolean edgeIsBackwards = lastVert > currentVert;
                numContactsDrawn++;
                ig2.setPaint(Color.BLACK);
                iSpatIndex = i;
                int j = i - 1;
                jSpatIndex = j;
                if (iSpatIndex < jSpatIndex) {
                    leftVertSpat = iSpatIndex;
                    leftVertSeq = currentVert;
                    rightVertSpat = jSpatIndex;
                    rightVertSeq = lastVert;
                } else {
                    leftVertSpat = jSpatIndex;
                    leftVertSeq = lastVert;
                    rightVertSpat = iSpatIndex;
                    rightVertSeq = currentVert;
                }
                leftVertPosX = vertStartX + (leftVertSeq * vertDist);
                rightVertPosX = vertStartX + (rightVertSeq * vertDist);
                connWidth = rightVertPosX - leftVertPosX;
                connHeight = connWidth / 2;
                connCenterX = rightVertPosX - (connWidth / 2);
                connCenterY = vertStartY - (connHeight / 2);
                connTopLeftX = leftVertPosX;
                connTopLeftY = vertStartY - (connHeight / 2);
                spacerX = vertWidth;
                spacerY = 0;
                if (Objects.equals(orientationsSeqOrder[leftVertSeq], FoldingGraph.ORIENTATION_UPWARDS)) {
                    leftVertPosY = vertStartY - vertHeight;
                    startUpwards = true;
                    if (debugDrawingKEY) {
                        DP.getInstance().d("  leftVert " + leftVertSeq + " starts upwards, ");
                    }
                } else {
                    leftVertPosY = vertStartY;
                    startUpwards = false;
                    if (debugDrawingKEY) {
                        DP.getInstance().d("  leftVert " + leftVertSeq + " starts downwards, ");
                    }
                }
                if (Objects.equals(orientationsSeqOrder[rightVertSeq], FoldingGraph.ORIENTATION_UPWARDS)) {
                    rightVertPosY = vertStartY;
                    if (debugDrawingKEY) {
                        DP.getInstance().d("  rightVert " + rightVertSeq + " starts upwards. ");
                    }
                } else {
                    rightVertPosY = vertStartY - vertHeight;
                    if (debugDrawingKEY) {
                        DP.getInstance().d("  rightVert " + rightVertSeq + " starts downwards. ");
                    }
                }
                boolean edgeIsCrossOver = Objects.equals(orientationsSeqOrder[leftVertSeq], orientationsSeqOrder[rightVertSeq]);
                if (debugDrawingKEY) {
                    System.out.print("Getting arc from " + leftVertPosX + "," + leftVertPosY + " to " + rightVertPosX + "," + rightVertPosY + ".\n");
                    Color old = ig2.getColor();
                    ig2.setColor(Color.BLACK);
                    int spacerDebugX = 10;
                    int spacerDebugY = i * 10;
                    ig2.drawString("i=" + i + "(" + lastVert + "->" + currentVert + ")", leftVertPosX + spacerDebugX, leftVertPosY + spacerDebugY);
                    ig2.drawString("x" + i, leftVertPosX, leftVertPosY);
                    ig2.drawString("y" + i, rightVertPosX, rightVertPosY);
                    ig2.setColor(old);
                }
                ArrayList<Shape> connShapes;
                /** If the central line of the crossover-connector would cut through a vertex, we may want to shift it to the side a bit. */
                int pixelsToShiftOnYAxis = 0;
                if (edgeIsCrossOver) {
                    if (relDrawDistToLast % 2 == 0) {
                        pixelsToShiftOnYAxis = (pl.vertDist / 4) + (relDrawDistToLast * 5);
                        if (pixelsToShiftOnYAxis > (pl.vertDist * 0.8)) {
                            pixelsToShiftOnYAxis = (pl.vertDist / 2);
                        }
                    }
                }
                if (edgeIsBackwards && edgeIsCrossOver) {
                    if (Settings.getBoolean("PTGLgraphComputation_B_key_use_alternate_arcs")) {
                        connShapes = DrawTools.getArcConnectorAlternative(leftVertPosX, rightVertPosY, rightVertPosX, leftVertPosY, ig2.getStroke(), !startUpwards, pixelsToShiftOnYAxis);
                    } else {
                        connShapes = DrawTools.getArcConnector(leftVertPosX, rightVertPosY, rightVertPosX, leftVertPosY, ig2.getStroke(), !startUpwards, pixelsToShiftOnYAxis);
                    }
                } else {
                    if (Settings.getBoolean("PTGLgraphComputation_B_key_use_alternate_arcs")) {
                        connShapes = DrawTools.getArcConnectorAlternative(leftVertPosX, leftVertPosY, rightVertPosX, rightVertPosY, ig2.getStroke(), startUpwards, pixelsToShiftOnYAxis);
                    } else {
                        connShapes = DrawTools.getArcConnector(leftVertPosX, leftVertPosY, rightVertPosX, rightVertPosY, ig2.getStroke(), startUpwards, pixelsToShiftOnYAxis);
                    }
                }
                for (Shape s : connShapes) {
                    ig2.draw(s);
                }
            }
            if (numContactsDrawn != fg.getEdgeList().size()) {
                if (fg.isASingleCycle() && numContactsDrawn == (fg.getEdgeList().size() - 1)) {
                } else {
                    DP.getInstance().w("SSEGraph", "drawFoldingGrapgKEYG2D: drew " + numContactsDrawn + ", but FG " + fg.toShortString() + " has only " + fg.getEdgeList().size() + " edges. " + (fg.isASingleCycle() ? " FG is single circle." : "FG is NOT single circle."));
                }
            }
        }
        ig2.setStroke(new BasicStroke(1));
        Polygon p;
        AffineTransform origXform = ig2.getTransform();
        AffineTransform newXform;
        int rotationCenterX;
        int rotationCenterY;
        int angle = 180;
        Shape shape;
        for (Integer s = 0; s < fg.getSize(); s++) {
            Integer i = fg.getSpatOrder().get(s);
            if (fg.getVertex(i).isHelix()) {
                ig2.setPaint(Color.RED);
            } else if (fg.getVertex(i).isBetaStrand()) {
                ig2.setPaint(Color.BLACK);
            } else if (fg.getVertex(i).isLigandSSE()) {
                ig2.setPaint(Color.MAGENTA);
            } else if (fg.getVertex(i).isOtherSSE()) {
                ig2.setPaint(Color.GRAY);
            } else {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            Integer currentVertX = vertStartX + (i * vertDist);
            Integer currentVertY = vertStartY;
            if (fg.getVertex(i).isHelix()) {
                p = DrawTools.getDefaultArrowPolygonUp(vertStartY - vertHeight, currentVertX, currentVertY);
            } else {
                p = DrawTools.getDefaultArrowPolygonUp(vertStartY - vertHeight, currentVertX, currentVertY);
            }
            shape = ig2.getStroke().createStrokedShape(p);
            newXform = (AffineTransform) (origXform.clone());
            rotationCenterX = (vertStartX + (i * vertDist));
            rotationCenterY = vertStartY - (vertHeight / 2);
            newXform.rotate(Math.toRadians(angle), rotationCenterX, rotationCenterY);
            if (Objects.equals(orientationsSeqOrder[i], FoldingGraph.ORIENTATION_DOWNWARDS)) {
                ig2.setTransform(newXform);
            }
            ig2.draw(shape);
            if (Objects.equals(orientationsSeqOrder[i], FoldingGraph.ORIENTATION_DOWNWARDS)) {
                ig2.setTransform(origXform);
            }
            Integer compareToTerminus = i;
            Integer ssePos = s;
            if (Objects.equals(s, 0)) {
                ig2.setPaint(Color.BLACK);
                ig2.drawString("N", currentVertX, currentVertY + 20);
            }
            if (Objects.equals(s, fg.numVertices() - 1)) {
                ig2.setPaint(Color.BLACK);
                ig2.drawString("C", currentVertX, currentVertY + 20);
            }
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {
            Integer printNth = 1;
            if (fg.getSize() > 0) {
                ig2.setPaint(Color.BLACK);
                ig2.drawString("FG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                ig2.setPaint(Color.LIGHT_GRAY);
                ig2.drawString("SQ", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                ig2.drawString("PG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + lineHeight + (stringHeight / 4));
                ig2.setPaint(Color.BLACK);
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            String sseNumberFoldingGraphKEY;
            for (Integer i = leftMostVertexInParent; i <= rightMostVertexInParent; i++) {
                Integer parentVertexPosInFoldingGraph = parentVertexPosInFG.get(i);
                if (parentVertexPosInFoldingGraph < 0) {
                    continue;
                }
                if ((i + 1) % printNth == 0) {
                    sseNumberFoldingGraph = "" + (parentVertexPosInFG.get(i) >= 0 ? (parentVertexPosInFG.get(i) + 1) : "");
                    sseNumberFoldingGraphKEY = fg.getSpatOrder().indexOf(i) + "";
                    sseNumberSeq = "" + (pg.getVertex(i).getSSESeqChainNum());
                    sseNumberProteinGraph = "" + (i + 1);
                    stringHeight = fontMetrics.getAscent();
                    ig2.setPaint(Color.BLACK);
                    ig2.drawString(sseNumberFoldingGraph, pl.getFooterStart().x + ((i - shiftBack[i]) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (stringHeight / 4));
                    ig2.setPaint(Color.LIGHT_GRAY);
                    ig2.drawString(sseNumberSeq, pl.getFooterStart().x + ((i - shiftBack[i]) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                    ig2.drawString(sseNumberProteinGraph, pl.getFooterStart().x + ((i - shiftBack[i]) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));
                    ig2.setPaint(Color.BLACK);
                }
            }
            if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
                ProteinGraphDrawer.drawLegendKEY(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4)), pl, fg);
            }
        }
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());
        DrawResult drawRes = new DrawResult(ig2, roi);
        return drawRes;
    }

    /**
     * Draws the RED folding graph image of this graph and returns the DrawResult.
     * If 'nonProteinGraph' is true, this graph is considered a custom (=non-protein) graph and the color coding for vertices and edges is NOT used.
     * In that case, the graph is drawn black and white and the labels for the N- and C-termini are NOT drawn.
     *
     * @param pnfr a folding graph notation result
     * @return the DrawResult. You can write this to a file or whatever.
     */
    private static DrawResult drawFoldingGraphREDG2D(PTGLNotationFoldResult pnfr) {
        FoldingGraph fg = pnfr.getFoldingGraph();
        SSEGraph pg = fg.getParent();
        Integer startVertexInParent = fg.getMinimalVertexIndexInParentGraph();
        Integer endVertexInParent = fg.getMaximalVertexIndexInParentGraph();
        List<Integer> fgVertexPosInParent = fg.getVertexIndexListInParentGraph();
        List<Integer> parentVertexPosInFG = ProtGraph.computeParentGraphVertexPositionsInFoldingGraph(fgVertexPosInParent, pg.getSize());
        Integer[] shiftBack = new Integer[pg.getSize()];
        Integer shift = 0;
        for (int i = 0; i < parentVertexPosInFG.size(); i++) {
            if (parentVertexPosInFG.get(i) < 0) {
                shift++;
            }
            shiftBack[i] = shift;
        }
        Integer numVerts = fg.getSize();
        Boolean bw = false;
        PageLayout pl = new PageLayout(numVerts);
        Position2D vertStart = pl.getVertStart();
        Integer lineHeight = pl.textLineHeight;
        SVGGraphics2D ig2;
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        ig2 = new SVGGraphics2D(document);
        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);
        Font font = pl.getStandardFont();
        Font fontBold = pl.getStandardFontBold();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();
        
        String proteinHeader = "RED " + pg.getGraphType() + " folding graph " + fg.getFoldingGraphFoldName() + " (# " + fg.getFoldingGraphNumber() + ") of PDB entry " + pg.getPdbid() + ", chain " + pg.getChainid() + "";
        proteinHeader += "\n" + Settings.getImagesCreatorVersionLine();
        if (FileParser.getWhichSseInfo().equals("author")){
            proteinHeader += "\n" + AUTHOR_SSE_WARNING;
        }
        
        String notation = "RED notation: '" + pnfr.redNotation + "'";
        Integer stringHeight = fontMetrics.getAscent();
        String sseNumberSeq;
        String sseNumberFoldingGraph;
        String sseNumberProteinGraph;
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            DrawTools.drawStringLineBreaks(ig2, proteinHeader, pl.headerStart.x, pl.headerStart.y);
            if (Settings.getBoolean("PTGLgraphComputation_B_print_notations_on_fg_images")) {
                ig2.drawString(notation, pl.headerStart.x, pl.headerStart.y + lineHeight);
            }
        }
        Integer k;
        Integer l;
        Shape shape;
        Arc2D.Double arc;
        ig2.setStroke(new BasicStroke(2));
        Integer edgeType;
        Integer leftVert;
        Integer rightVert;
        Integer leftVertPosX;
        Integer rightVertPosX;
        Integer arcWidth;
        Integer arcHeight;
        Integer arcTopLeftX;
        Integer arcTopLeftY;
        Integer spacerX;
        Integer spacerY;
        for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
            for (Integer j = i + 1; j <= endVertexInParent; j++) {
                k = parentVertexPosInFG.get(i);
                l = parentVertexPosInFG.get(j);
                if (k < 0 || l < 0) {
                    continue;
                }
                if (fg.containsEdge(k, l)) {
                    edgeType = fg.getContactSpatRel(k, l);
                    if (edgeType.equals(SpatRel.PARALLEL)) {
                        ig2.setPaint(Color.RED);
                    } else if (edgeType.equals(SpatRel.ANTIPARALLEL)) {
                        ig2.setPaint(Color.BLUE);
                    } else if (edgeType.equals(SpatRel.MIXED)) {
                        ig2.setPaint(Color.GREEN);
                    } else if (edgeType.equals(SpatRel.LIGAND)) {
                        ig2.setPaint(Color.MAGENTA);
                    } else if (edgeType.equals(SpatRel.BACKBONE)) {
                        ig2.setPaint(Color.ORANGE);
                    } else {
                        ig2.setPaint(Color.LIGHT_GRAY);
                    }
                    if (bw) {
                        ig2.setPaint(Color.LIGHT_GRAY);
                    }
                    if (i < j) {
                        leftVert = i;
                        rightVert = j;
                    } else {
                        leftVert = j;
                        rightVert = i;
                    }
                    leftVertPosX = pl.getVertStart().x + ((leftVert - shiftBack[leftVert]) * pl.vertDist);
                    rightVertPosX = pl.getVertStart().x + ((rightVert - shiftBack[rightVert]) * pl.vertDist);
                    arcWidth = rightVertPosX - leftVertPosX;
                    arcHeight = arcWidth / 2;
                    arcTopLeftX = leftVertPosX;
                    arcTopLeftY = pl.getVertStart().y - arcHeight / 2;
                    spacerX = pl.vertRadius;
                    spacerY = 0;
                    arc = new Arc2D.Double(arcTopLeftX + spacerX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
                    shape = ig2.getStroke().createStrokedShape(arc);
                    ig2.fill(shape);
                }
            }
        }
        Ellipse2D.Double circle;
        Rectangle2D.Double rect;
        ig2.setStroke(new BasicStroke(2));
        for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
            if (pg.getVertex(i).isHelix()) {
                ig2.setPaint(Color.RED);
            } else if (pg.getVertex(i).isBetaStrand()) {
                ig2.setPaint(Color.BLACK);
            } else if (pg.getVertex(i).isLigandSSE()) {
                ig2.setPaint(Color.MAGENTA);
            } else if (pg.getVertex(i).isOtherSSE()) {
                ig2.setPaint(Color.GRAY);
            } else {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            Integer parentVertexPosInFoldingGraph = parentVertexPosInFG.get(i);
            if (parentVertexPosInFoldingGraph < 0) {
                continue;
            }
            if (bw) {
                ig2.setPaint(Color.GRAY);
            }
            SSE sse = pg.getVertex(i);
            if (sse.isBetaStrand()) {
                rect = new Rectangle2D.Double(vertStart.x + ((i - shiftBack[i]) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(rect);
            } else if (sse.isLigandSSE()) {
                circle = new Ellipse2D.Double(vertStart.x + ((i - shiftBack[i]) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.setStroke(new BasicStroke(3));
                ig2.draw(circle);
                ig2.setStroke(new BasicStroke(2));
            } else {
                circle = new Ellipse2D.Double(vertStart.x + ((i - shiftBack[i]) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(circle);
            }
        }
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        if (!bw) {
            if (fg.getSize() > 0) {
                ig2.setFont(fontBold);
                ig2.drawString("N", vertStart.x - pl.vertDist, vertStart.y + 20);
                ig2.drawString("C", vertStart.x + numVerts * pl.vertDist, vertStart.y + 20);
                ig2.setFont(font);
            }
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {
            Integer printNth = 1;
            if (fg.getSize() > 99) {
                printNth = 2;
            }
            if (fg.getSize() > 999) {
                printNth = 3;
            }
            if (fg.getSize() > 0) {
                ig2.drawString("FG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                ig2.setPaint(Color.LIGHT_GRAY);
                ig2.drawString("SQ", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                ig2.drawString("PG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + lineHeight + (stringHeight / 4));
                ig2.setPaint(Color.BLACK);
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
                Integer parentVertexPosInFoldingGraph = parentVertexPosInFG.get(i);
                if (parentVertexPosInFoldingGraph < 0) {
                    continue;
                }
                if ((i + 1) % printNth == 0) {
                    sseNumberFoldingGraph = "" + (parentVertexPosInFG.get(i) >= 0 ? (parentVertexPosInFG.get(i) + 1) : "");
                    sseNumberSeq = "" + (pg.getVertex(i).getSSESeqChainNum());
                    sseNumberProteinGraph = "" + (i + 1);
                    stringHeight = fontMetrics.getAscent();
                    ig2.setPaint(Color.BLACK);
                    ig2.drawString(sseNumberFoldingGraph, pl.getFooterStart().x + ((i - shiftBack[i]) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (stringHeight / 4));
                    ig2.setPaint(Color.LIGHT_GRAY);
                    ig2.drawString(sseNumberSeq, pl.getFooterStart().x + ((i - shiftBack[i]) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                    ig2.drawString(sseNumberProteinGraph, pl.getFooterStart().x + ((i - shiftBack[i]) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));
                    ig2.setPaint(Color.BLACK);
                }
            }
            if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
                ProteinGraphDrawer.drawLegend(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4)), pl, fg);
            }
        }
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());
        DrawResult drawRes = new DrawResult(ig2, roi);
        return drawRes;
    }

    /**
     * Draws the DEF folding graph image of this graph and returns the DrawResult.
     *
     * The DEF folding graph visualization was introduced by Tim. It is similar to ADJ, but draws all vertices of the parent PG; not just those between the start and end vertices of the FG.
     * Note that this does not have a linnot string, it is just for visualization of the FG within the PG.
     *
     * If 'nonProteinGraph' is true, this graph is considered a custom (=non-protein) graph and the color coding for vertices and edges is NOT used.
     * In that case, the graph is drawn black and white and the labels for the N- and C-termini are NOT drawn.
     *
     * @param pnfr a folding graph notation result
     * @return the DrawResult. You can write this to a file or whatever.
     */
    private static DrawResult drawFoldingGraphDEFG2D(PTGLNotationFoldResult pnfr) {
        FoldingGraph fg = pnfr.getFoldingGraph();
        SSEGraph pg = fg.getParent();
        Integer startVertexInParent = 0;
        Integer endVertexInParent = pg.getSize() - 1;
        Integer shiftBack = startVertexInParent;
        List<Integer> fgVertexPosInParent = fg.getVertexIndexListInParentGraph();
        List<Integer> parentVertexPosInFG = pg.computeParentGraphVertexPositionsInFoldingGraph(fgVertexPosInParent, pg.getSize());
        Integer numVerts = endVertexInParent - startVertexInParent + 1;
        Boolean bw = false;
        PageLayout pl = new PageLayout(numVerts);
        Position2D vertStart = pl.getVertStart();
        Integer lineHeight = pl.textLineHeight;
        SVGGraphics2D ig2;
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        ig2 = new SVGGraphics2D(document);
        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);
        Font font = pl.getStandardFont();
        Font fontBold = pl.getStandardFontBold();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();
        
        String proteinHeader = "DEF " + pg.getGraphType() + " folding graph " + fg.getFoldingGraphFoldName() + " (# " + fg.getFoldingGraphNumber() + ") of PDB entry " + pg.getPdbid() + ", chain " + pg.getChainid() + "";
        proteinHeader += "\n" + Settings.getImagesCreatorVersionLine();
        if (FileParser.getWhichSseInfo().equals("author")){
            proteinHeader += "\n" + AUTHOR_SSE_WARNING;
        }
        
        Integer stringHeight = fontMetrics.getAscent();
        String sseNumberSeq;
        String sseNumberFoldingGraph;
        String sseNumberProteinGraph;
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            DrawTools.drawStringLineBreaks(ig2, proteinHeader, pl.headerStart.x, pl.headerStart.y);
        }
        Integer k;
        Integer l;
        Shape shape;
        Arc2D.Double arc;
        ig2.setStroke(new BasicStroke(2));
        Integer edgeType;
        Integer leftVert;
        Integer rightVert;
        Integer leftVertPosX;
        Integer rightVertPosX;
        Integer arcWidth;
        Integer arcHeight;
        Integer arcTopLeftX;
        Integer arcTopLeftY;
        Integer spacerX;
        Integer spacerY;
        for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
            for (Integer j = i + 1; j <= endVertexInParent; j++) {
                Boolean parentGraphOnlyPair;
                k = parentVertexPosInFG.get(i);
                l = parentVertexPosInFG.get(j);
                if (k < 0 || l < 0) {
                    parentGraphOnlyPair = true;
                } else {
                    parentGraphOnlyPair = false;
                }
                if (parentGraphOnlyPair) {
                    if (pg.containsEdge(i, j)) {
                        ig2.setPaint(Color.LIGHT_GRAY);
                        if (i < j) {
                            leftVert = i;
                            rightVert = j;
                        } else {
                            leftVert = j;
                            rightVert = i;
                        }
                        leftVertPosX = pl.getVertStart().x + ((leftVert - shiftBack) * pl.vertDist);
                        rightVertPosX = pl.getVertStart().x + ((rightVert - shiftBack) * pl.vertDist);
                        arcWidth = rightVertPosX - leftVertPosX;
                        arcHeight = arcWidth / 2;
                        arcTopLeftX = leftVertPosX;
                        arcTopLeftY = pl.getVertStart().y - arcHeight / 2;
                        spacerX = pl.vertRadius;
                        spacerY = 0;
                        arc = new Arc2D.Double(arcTopLeftX + spacerX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
                        shape = ig2.getStroke().createStrokedShape(arc);
                        ig2.fill(shape);
                    }
                } else {
                    if (fg.containsEdge(k, l)) {
                        edgeType = fg.getContactSpatRel(k, l);
                        if (edgeType.equals(SpatRel.PARALLEL)) {
                            ig2.setPaint(Color.RED);
                        } else if (edgeType.equals(SpatRel.ANTIPARALLEL)) {
                            ig2.setPaint(Color.BLUE);
                        } else if (edgeType.equals(SpatRel.MIXED)) {
                            ig2.setPaint(Color.GREEN);
                        } else if (edgeType.equals(SpatRel.LIGAND)) {
                            ig2.setPaint(Color.MAGENTA);
                        } else if (edgeType.equals(SpatRel.BACKBONE)) {
                            ig2.setPaint(Color.ORANGE);
                        } else {
                            ig2.setPaint(Color.LIGHT_GRAY);
                        }
                        if (bw) {
                            ig2.setPaint(Color.LIGHT_GRAY);
                        }
                        if (i < j) {
                            leftVert = i;
                            rightVert = j;
                        } else {
                            leftVert = j;
                            rightVert = i;
                        }
                        leftVertPosX = pl.getVertStart().x + ((leftVert - shiftBack) * pl.vertDist);
                        rightVertPosX = pl.getVertStart().x + ((rightVert - shiftBack) * pl.vertDist);
                        arcWidth = rightVertPosX - leftVertPosX;
                        arcHeight = arcWidth / 2;
                        arcTopLeftX = leftVertPosX;
                        arcTopLeftY = pl.getVertStart().y - arcHeight / 2;
                        spacerX = pl.vertRadius;
                        spacerY = 0;
                        arc = new Arc2D.Double(arcTopLeftX + spacerX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
                        shape = ig2.getStroke().createStrokedShape(arc);
                        ig2.fill(shape);
                    }
                }
            }
        }
        Ellipse2D.Double circle;
        Rectangle2D.Double rect;
        ig2.setStroke(new BasicStroke(2));
        for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
            if (pg.getVertex(i).isHelix()) {
                ig2.setPaint(Color.RED);
            } else if (pg.getVertex(i).isBetaStrand()) {
                ig2.setPaint(Color.BLACK);
            } else if (pg.getVertex(i).isLigandSSE()) {
                ig2.setPaint(Color.MAGENTA);
            } else if (pg.getVertex(i).isOtherSSE()) {
                ig2.setPaint(Color.GRAY);
            } else {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            Integer parentVertexPosInFoldingGraph = parentVertexPosInFG.get(i);
            if (parentVertexPosInFoldingGraph < 0) {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            if (bw) {
                ig2.setPaint(Color.GRAY);
            }
            SSE sse = pg.getVertex(i);
            if (sse.isBetaStrand()) {
                rect = new Rectangle2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(rect);
            } else if (sse.isLigandSSE()) {
                circle = new Ellipse2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.setStroke(new BasicStroke(3));
                ig2.draw(circle);
                ig2.setStroke(new BasicStroke(2));
            } else {
                circle = new Ellipse2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(circle);
            }
        }
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        if (!bw) {
            if (fg.getSize() > 0) {
                ig2.setFont(fontBold);
                ig2.drawString("N", vertStart.x - pl.vertDist, vertStart.y + 20);
                ig2.drawString("C", vertStart.x + numVerts * pl.vertDist, vertStart.y + 20);
                ig2.setFont(font);
            }
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {
            Integer printNth = 1;
            if (fg.getSize() > 99) {
                printNth = 2;
            }
            if (fg.getSize() > 999) {
                printNth = 3;
            }
            if (fg.getSize() > 0) {
                ig2.drawString("FG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                ig2.setPaint(Color.LIGHT_GRAY);
                ig2.drawString("SQ", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                ig2.drawString("PG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + lineHeight + (stringHeight / 4));
                ig2.setPaint(Color.BLACK);
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
                if ((i + 1) % printNth == 0) {
                    sseNumberFoldingGraph = "" + (parentVertexPosInFG.get(i) >= 0 ? (parentVertexPosInFG.get(i) + 1) : "");
                    sseNumberSeq = "" + (pg.getVertex(i).getSSESeqChainNum());
                    sseNumberProteinGraph = "" + (i + 1);
                    stringHeight = fontMetrics.getAscent();
                    ig2.setPaint(Color.BLACK);
                    ig2.drawString(sseNumberFoldingGraph, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (stringHeight / 4));
                    ig2.setPaint(Color.LIGHT_GRAY);
                    ig2.drawString(sseNumberSeq, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                    ig2.drawString(sseNumberProteinGraph, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));
                    ig2.setPaint(Color.BLACK);
                }
            }
            if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
                ProteinGraphDrawer.drawLegendDEFADJ(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4)), pl, fg, true);
            }
        }
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());
        DrawResult drawRes = new DrawResult(ig2, roi);
        return drawRes;
    }

    /**
     * Draws the legend for the graph at the given position. This legend is not suitable for SEQ folding graphs, because their edges are different.
     * This is a special version for DEF and ADJ graphs, which adds the gray vertices and labels them as parent graph vertices.
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position (x, y) where to start drawing
     * @param pl the page layout
     * @param g tge sse graph
     * @param includeParentEdge whether to draw the symbol fro parent edges (only DEF graphs need this, but not ADJ graphs)
     * @return the x coordinate in the image where the legend ends (which is the left margin + the legend width).
     * This can be used to determine the minimal width of the total image (it has to be at least this value).
     */
    public static Integer drawLegendDEFADJ(SVGGraphics2D ig2, Position2D startPos, PageLayout pl, SSEGraph g, Boolean includeParentEdge) {
        Boolean drawAll = Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend_always_all");
        ig2.setFont(pl.getLegendFont());
        FontMetrics fontMetrics = ig2.getFontMetrics();
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        Integer spacer = 10;
        Integer pixelPosX = startPos.x;
        Integer vertWidth = pl.getVertDiameter();
        Integer vertOffset = pl.getVertDiameter() / 4 * 3;
        String label;
        if (g.numEdges() > 0) {
            label = "[Edges: ";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsContactTypeParallel() || drawAll) {
            label = "parallel";
            ig2.setPaint(Color.RED);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsContactTypeAntiparallel() || drawAll) {
            label = "antiparallel";
            ig2.setPaint(Color.BLUE);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsContactTypeMixed() || drawAll) {
            label = "mixed";
            ig2.setPaint(Color.GREEN);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsContactTypeLigand() || drawAll) {
            label = "ligand";
            ig2.setPaint(Color.MAGENTA);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        label = "parent";
        ig2.setPaint(Color.GRAY);
        ig2.drawString(label, pixelPosX, startPos.y);
        pixelPosX += fontMetrics.stringWidth(label) + spacer;
        if (g.getChainEnds().size() > 0 || drawAll) {
            label = "interchain";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.numEdges() > 0) {
            label = "]";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.numVertices() > 0) {
            label = " [Vertices: ";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeHelix() || drawAll) {
            ProteinGraphDrawer.drawSymbolAlphaHelix(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "helix";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeBetaStrand() || drawAll) {
            ProteinGraphDrawer.drawSymbolBetaStrand(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "strand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeLigand() || drawAll) {
            ProteinGraphDrawer.drawSymbolLigand(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "ligand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeOther() || drawAll) {
            ProteinGraphDrawer.drawSymbolOtherSSE(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "other";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        ProteinGraphDrawer.drawSymbolsParentSSEs(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
        pixelPosX += vertWidth + spacer + 50;
        label = "parent";
        ig2.drawString(label, pixelPosX, startPos.y);
        pixelPosX += fontMetrics.stringWidth(label) + spacer;
        if (g.numVertices() > 0) {
            label = "]";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        return pixelPosX;
    }

    
    /**
     * Draws the protein graph image of this graph and returns the DrawResult.
     * If 'nonProteinGraph' is true, this graph is considered a custom (=non-protein) graph and the color coding for vertices and edges is NOT used.
     * In that case, the graph is drawn black and white and the labels for the N- and C-termini are NOT drawn.
     *
     * @param nonProteinGraph whether the graph is a non-protein graph and thus does NOT contain information on the relative SSE orientation in the expected way. If so, it is drawn in gray scale because the color code becomes useless (true => gray scale, false => color).
     * @param pg the graph to draw
     * @return the DrawResult. You can write this to a file or whatever.
     */
    private static DrawResult drawProteinGraphG2D(Boolean nonProteinGraph, ProtGraph pg) {
        Map<Integer, String> vertexMarkings = new HashMap<>();
        List<String> ignoreChains = new ArrayList<>();
        return ProteinGraphDrawer.drawProteinGraphG2D(nonProteinGraph, pg, vertexMarkings, ignoreChains);
    }
    
    /**
     * Draws the protein graph image of this graph and returns the DrawResult.
     * If 'nonProteinGraph' is true, this graph is considered a custom (=non-protein) graph and the color coding for vertices and edges is NOT used.
     * In that case, the graph is drawn black and white and the labels for the N- and C-termini are NOT drawn.
     *
     * @param nonProteinGraph whether the graph is a non-protein graph and thus does NOT contain information on the relative SSE orientation in the expected way. If so, it is drawn in gray scale because the color code becomes useless (true => gray scale, false => color).
     * @param pg the graph to draw
     * @param vertexMarkings a map of special markings for vertices, supply an empty Map if no vertices should be marked in output image. This can be used to visually emphasize subsets of vertices in the graph.
     * @return the DrawResult. You can write this to a file or whatever.
     */
    private static DrawResult drawProteinGraphG2D(Boolean nonProteinGraph, ProtGraph pg, Map<Integer, String> vertexMarkings, List<String> ignoreChains) {
        
        // generate a list of ignored SSEs from the list of ignored chains
        List<Integer> ignoredSSEIndices = new ArrayList<>();
        for(int i = 0; i < pg.getSize(); i++) {
            SSE s = pg.getVertex(i);
            if(ignoreChains.contains(pg.getChainNameOfSSE(i))) {
                ignoredSSEIndices.add(i);
            }
        }
        int numDrawnVertices = pg.getSize() - ignoredSSEIndices.size();
                
        Integer numVerts = pg.numVertices();
        Boolean bw = nonProteinGraph;
        PageLayout pl = new PageLayout(numDrawnVertices);
        Position2D vertStart = pl.getVertStart();
        SVGGraphics2D ig2;
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        ig2 = new SVGGraphics2D(document);
        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);
        Font font = pl.getStandardFont();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();
        
        String proteinHeader = "" + pg.getGraphType() + " Protein Graph of PDB entry " + pg.getPdbid() + ", chain " + pg.getChainid() + "";
        proteinHeader += "\n" + Settings.getImagesCreatorVersionLine();
        if (FileParser.getWhichSseInfo().equals("author")){
            proteinHeader += "\n" + AUTHOR_SSE_WARNING;
        }
        
        Integer stringHeight = fontMetrics.getAscent();
        String sseNumberSeq;
        String sseNumberGraph;
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            DrawTools.drawStringLineBreaks(ig2, proteinHeader, pl.headerStart.x, pl.headerStart.y);
        }
        Shape shape;
        Arc2D.Double arc;
        ig2.setStroke(new BasicStroke(2));
        Integer edgeType;
        Integer leftVert;
        Integer rightVert;
        Integer leftVertPosX;
        Integer rightVertPosX;
        Integer arcWidth;
        Integer arcHeight;
        Integer arcTopLeftX;
        Integer arcTopLeftY;
        Integer spacerX;
        Integer spacerY;
        Integer iChainID;
        Integer jChainID;
        
        Color C_PARALLEL = Color.RED;
        Color C_MIXED = Color.GREEN;
        Color C_ANTIPARALLEL = Color.BLUE;
        Color C_LIGAND = Color.MAGENTA;
        Color C_BACKBONE = Color.ORANGE;
        Color C_OTHER = Color.LIGHT_GRAY;
        
        Color C_INTERCHAIN_PARALLEL = PlccUtilities.mutateColor(C_PARALLEL, 50, 50, 50);
        Color C_INTERCHAIN_MIXED = PlccUtilities.mutateColor(C_MIXED, 50, 50, 50);
        Color C_INTERCHAIN_ANTIPARALLEL = PlccUtilities.mutateColor(C_ANTIPARALLEL, 50, 50, 50);
        Color C_INTERCHAIN_LIGAND = PlccUtilities.mutateColor(C_LIGAND, 50, 50, 50);
        Color C_INTERCHAIN_BACKBONE = PlccUtilities.mutateColor(C_BACKBONE, 50, 50, 50); // contact type does not really make sense (backbone interchain)
        Color C_INTERCHAIN_OTHER = PlccUtilities.mutateColor(C_OTHER, 50, 50, 50);
        
        Color C_INTERCHAIN_DIFFMOL_PARALLEL = PlccUtilities.mutateColor(C_PARALLEL, 100, 100, 100);
        Color C_INTERCHAIN_DIFFMOL_MIXED = PlccUtilities.mutateColor(C_MIXED, 100, 100, 100);
        Color C_INTERCHAIN_DIFFMOL_ANTIPARALLEL = PlccUtilities.mutateColor(C_ANTIPARALLEL, 100, 100, 100);
        Color C_INTERCHAIN_DIFFMOL_LIGAND = PlccUtilities.mutateColor(C_LIGAND, 100, 100, 100);
        Color C_INTERCHAIN_DIFFMOL_BACKBONE = PlccUtilities.mutateColor(C_BACKBONE, 100, 100, 100); // contact type does not really make sense (backbone interchain)
        Color C_INTERCHAIN_DIFFMOL_OTHER = PlccUtilities.mutateColor(C_OTHER, 100, 100, 100);
        
        
        
        //DP.getInstance().i("ProteinGraphDrawer", "Ignoring " + ignoreChains.size() + " chains (" + ignoredSSEIndices.size() + " of " + pg.getSize() + " SSEs).");
        
        // Compute the draw position array. This defines the position of each vertex (given by index in the graph) in the drawing order (the drawing order may be  shifted to the left due to ignored vertices, so that SSE n is drawn at a position < n).
        int[] drawPos = new int[pg.getSize()];
        int shiftLeft = 0;
        final int DP_NONE = -1;
        for(int i = 0; i < pg.getSize(); i++) {
            if(ignoredSSEIndices.contains(i)) {
                drawPos[i] = DP_NONE;
                shiftLeft++;
            }
            else {                
                drawPos[i] = (i - shiftLeft);
            }
        }
        
        //DP.getInstance().i("ProteinGraphDrawer", "The drawPos array is: '" + IO.intArrayToString(drawPos) + "'.");
        
        float[] dashPattern = { 10, 20, 10, 20 }; // dash pattern for dashed lines: alternating lengths of transparent (even index) and opaque (uneven index) line segement lengths to draw
        float[] dashPatternDiffMol = { 20, 40, 20, 40 };
        
        for (Integer i = 0; i < pg.getSize(); i++) {
            for (Integer j = i + 1; j < pg.getSize(); j++) {
                if (pg.containsEdge(i, j)) {
                    
                    // skip edges of ignored chains
                    if(ignoreChains.contains(pg.getChainNameOfSSE(i)) || ignoreChains.contains(pg.getChainNameOfSSE(j))) {
                        continue;
                    }
                    
                    edgeType = pg.getContactSpatRel(i, j);
                    if (edgeType.equals(SpatRel.PARALLEL)) {
                        ig2.setPaint(C_PARALLEL);
                    } else if (edgeType.equals(SpatRel.ANTIPARALLEL)) {
                        ig2.setPaint(C_ANTIPARALLEL);
                    } else if (edgeType.equals(SpatRel.MIXED)) {
                        ig2.setPaint(C_MIXED);
                    } else if (edgeType.equals(SpatRel.LIGAND)) {
                        ig2.setPaint(C_LIGAND);
                    } else if (edgeType.equals(SpatRel.BACKBONE)) {
                        ig2.setPaint(C_BACKBONE);
                    } else if (edgeType.equals(SpatRel.COMPLEX)) {
                        ig2.setPaint(Color.BLACK);  // not really used, the color for interchain gets mutated from the spatRel color below            
                    } else {
                        ig2.setPaint(C_OTHER);
                    }
                    if (bw) {
                        ig2.setPaint(C_OTHER);
                    }
                                        
                    
                    iChainID = -1;
                    jChainID = -1;
                    String iChainMolID, jChainMolID;
                    for (Integer x = 0; x < pg.getChainEnds().size(); x++) {
                        if (i < pg.getChainEnds().get(x)) {
                            iChainID = x;
                            break;
                        }
                    }
                    
                    for (Integer x = 0; x < pg.getChainEnds().size(); x++) {
                        if (j < pg.getChainEnds().get(x)) {
                            jChainID = x;
                            break;
                        }
                    }
                    
                    // ------- set color for interchain (complex) contacts ------------
                    if (!Objects.equals(iChainID, jChainID)) {
                        
                        // interchain contact!
                        
                        // check whether this is an interchain contact between SSEs which belong to different macromolecules
                        iChainMolID = pg.getMolIDOfSSE(i);
                        jChainMolID = pg.getMolIDOfSSE(j);
                        if( ! Objects.equals(iChainMolID, jChainMolID)) {
                            //System.out.println("CGIC-D: Interchain contact between SSEs " + i + " and " + j + " is a contact between chains " + pg.getAllChains().get(iChainID).getPdbChainID() + " and " + pg.getAllChains().get(jChainID).getPdbChainID() + " of different macromolecules, MOL_IDs " + iChainMolID + " and " + jChainMolID + ".");
                            //ig2.setStroke(new BasicStroke(2, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 10, dashPattern, 0));
                            if (edgeType.equals(SpatRel.PARALLEL)) {
                                ig2.setPaint(C_INTERCHAIN_DIFFMOL_PARALLEL);
                            } else if (edgeType.equals(SpatRel.ANTIPARALLEL)) {
                                ig2.setPaint(C_INTERCHAIN_DIFFMOL_ANTIPARALLEL);
                            } else if (edgeType.equals(SpatRel.MIXED)) {
                                ig2.setPaint(C_INTERCHAIN_DIFFMOL_MIXED);
                            } else if (edgeType.equals(SpatRel.LIGAND)) {
                                ig2.setPaint(C_INTERCHAIN_DIFFMOL_LIGAND);
                            } else if (edgeType.equals(SpatRel.BACKBONE)) {
                                DP.getInstance().w("ProteinGraphDrawer", "Interchain backbone contact found, this makes no sense.");
                                ig2.setPaint(C_INTERCHAIN_DIFFMOL_BACKBONE);
                            } else {
                                ig2.setPaint(C_INTERCHAIN_DIFFMOL_OTHER);
                            }
                        }
                        else {
                            //System.out.println("CGIC-S: Interchain contact between SSEs " + i + " and " + j + " is a contact between chains " + pg.getAllChains().get(iChainID).getPdbChainID() + " and " + pg.getAllChains().get(jChainID).getPdbChainID() + " belonging to the same macromolecule (MOL_IDs " + iChainMolID + " and " + iChainMolID + ").");
                            //ig2.setStroke(new BasicStroke(2, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 10, dashPatternDiffMol, 0));
                            if (edgeType.equals(SpatRel.PARALLEL)) {
                                ig2.setPaint(C_INTERCHAIN_PARALLEL);
                            } else if (edgeType.equals(SpatRel.ANTIPARALLEL)) {
                                ig2.setPaint(C_INTERCHAIN_ANTIPARALLEL);
                            } else if (edgeType.equals(SpatRel.MIXED)) {
                                ig2.setPaint(C_INTERCHAIN_MIXED);
                            } else if (edgeType.equals(SpatRel.LIGAND)) {
                                ig2.setPaint(C_INTERCHAIN_LIGAND);
                            } else if (edgeType.equals(SpatRel.BACKBONE)) {
                                DP.getInstance().w("ProteinGraphDrawer", "Interchain backbone contact found, this makes no sense.");
                                ig2.setPaint(C_INTERCHAIN_BACKBONE);
                            } else {
                                ig2.setPaint(C_INTERCHAIN_OTHER);
                            }
                        }
                        
                        
                    }
                    else {
                        // intrachain contact
                        ig2.setStroke(new BasicStroke(2));
                    }
                    
                    
                    if (i < j) {
                        leftVert = i;
                        rightVert = j;
                    } else {
                        leftVert = j;
                        rightVert = i;
                    }
                    leftVertPosX = pl.getVertStart().x + (drawPos[leftVert] * pl.vertDist);
                    rightVertPosX = pl.getVertStart().x + (drawPos[rightVert] * pl.vertDist);
                    arcWidth = rightVertPosX - leftVertPosX;
                    arcHeight = arcWidth / 2;
                    arcTopLeftX = leftVertPosX;
                    arcTopLeftY = pl.getVertStart().y - arcHeight / 2;
                    spacerX = pl.vertRadius;
                    spacerY = 0;
                    arc = new Arc2D.Double(arcTopLeftX + spacerX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
                    shape = ig2.getStroke().createStrokedShape(arc);
                    ig2.draw(shape);
                }
            }
        }
        Ellipse2D.Double circle, markingCircle;
        Rectangle2D.Double rect, markingRect;
        ig2.setStroke(new BasicStroke(2));
        for (Integer i = 0; i < pg.getSize(); i++) {
            
            // skip vertices of ignored chains
            if(ignoreChains.contains(pg.getChainNameOfSSE(i))) {
                continue;
            }
            
            if (pg.getVertex(i).isHelix()) {
                ig2.setPaint(Color.RED);
            } else if (pg.getVertex(i).isBetaStrand()) {
                ig2.setPaint(Color.BLACK);
            } else if (pg.getVertex(i).isLigandSSE()) {
                ig2.setPaint(Color.MAGENTA);
            } else if (pg.getVertex(i).isOtherSSE()) {
                ig2.setPaint(Color.GRAY);
            } else {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            if (bw) {
                ig2.setPaint(Color.GRAY);
            }
            
            String vMarking = vertexMarkings.get(i);
            Color markingBorderColor = Color.GRAY;
            Color markingLabelColor = Color.GRAY;
            int markingWidth = 4;
            int markingOffsetY = 25;
            
            if (pg.getVertex(i).isBetaStrand()) {
                rect = new Rectangle2D.Double(vertStart.x + (drawPos[i] * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());                                
                ig2.fill(rect);
                if(vMarking != null) {
                    //System.out.println("#####Found marking for strand vertex " + i + ".");
                    markingRect = new Rectangle2D.Double(vertStart.x + (drawPos[i] * pl.vertDist), vertStart.y, pl.getVertDiameter() + (0 * markingWidth), pl.getVertDiameter() + (0 * markingWidth));
                    ig2.setColor(markingBorderColor);
                    ig2.setStroke(new BasicStroke(markingWidth));
                    ig2.draw(markingRect);
                    ig2.setColor(markingLabelColor);
                    ig2.drawString(vMarking, new Double(markingRect.getMinX()).intValue() + (pl.vertRadius / 2), new Double(markingRect.getCenterY()).intValue() + markingOffsetY);
                }
            } else if (pg.getVertex(i).isLigandSSE()) {
                int outerStrokeWidth = 4;
                int halfouterStrokeWidth = outerStrokeWidth / 2;
                circle = new Ellipse2D.Double(vertStart.x + (drawPos[i] * pl.vertDist) + halfouterStrokeWidth, vertStart.y + halfouterStrokeWidth, pl.getVertDiameter() - outerStrokeWidth, pl.getVertDiameter() - outerStrokeWidth);
                ig2.setStroke(new BasicStroke(outerStrokeWidth));
                ig2.draw(circle);
                ig2.setStroke(new BasicStroke(2));
                if(vMarking != null) {
                    //System.out.println("#####Found marking for ligand vertex " + i + ".");
                    markingCircle = new Ellipse2D.Double(vertStart.x + (drawPos[i] * pl.vertDist), vertStart.y, pl.getVertDiameter() + (0 * markingWidth), pl.getVertDiameter() + (0 * markingWidth));
                    ig2.setColor(markingBorderColor);
                    ig2.setStroke(new BasicStroke(markingWidth));
                    ig2.draw(markingCircle);
                    ig2.setColor(markingLabelColor);
                    ig2.drawString(vMarking, new Double(markingCircle.getMinX()).intValue() + (pl.vertRadius / 2), new Double(markingCircle.getCenterY()).intValue() + markingOffsetY);
                }
            } else {
                circle = new Ellipse2D.Double(vertStart.x + (drawPos[i] * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(circle);
                if(vMarking != null) {
                    //System.out.println("#####Found marking for other vertex " + i + ".");
                    markingCircle = new Ellipse2D.Double(vertStart.x + (drawPos[i] * pl.vertDist), vertStart.y, pl.getVertDiameter() + (0 * markingWidth), pl.getVertDiameter() + (0 * markingWidth));
                    ig2.setColor(markingBorderColor);
                    ig2.setStroke(new BasicStroke(markingWidth));
                    ig2.draw(markingCircle);
                    ig2.setColor(markingLabelColor);
                    ig2.drawString(vMarking, new Double(markingCircle.getMinX()).intValue() + (pl.vertRadius / 2), new Double(markingCircle.getCenterY()).intValue() + markingOffsetY);
                }
            }
        }
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        if (!bw) {
            if (pg.getSize() > 0 && (! pg.isComplexGraph())) {
                ig2.drawString("N", vertStart.x - pl.vertDist, vertStart.y + 20);
                ig2.drawString("C", vertStart.x + numDrawnVertices * pl.vertDist, vertStart.y + 20);
            }
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {
            Integer printNth = 1;
            if (pg.getSize() > 9) {
                printNth = 1;
            }
            if (pg.getSize() > 99) {
                printNth = 1;
            }
            if (pg.getSize() > 999) {
                printNth = 2;
            }
            Integer lineHeight = pl.textLineHeight;
            if (pg.getSize() > 0) {
                ig2.drawString("PG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                ig2.drawString("SQ", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            iChainID = -1;
            for (Integer i = 0; i < pg.getSize(); i++) {
                if ((i + 1) % printNth == 0) {
                    sseNumberGraph = "" + (drawPos[i] + 1);
                    sseNumberSeq = "" + (pg.getVertex(i).getSSESeqChainNum());
                    stringHeight = fontMetrics.getAscent();
                    if(drawPos[i] >= 0){
                        ig2.drawString(sseNumberGraph, pl.getFooterStart().x + (drawPos[i] * pl.vertDist) + (pl.vertRadius / 2), pl.getFooterStart().y + (stringHeight / 4));
                        ig2.drawString(sseNumberSeq, pl.getFooterStart().x + (drawPos[i] * pl.vertDist) + (pl.vertRadius / 2), pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                    }
                    for (Integer x = 0; x < pg.getChainEnds().size(); x++) {
                        if (i < pg.getChainEnds().get(x)) {
                            iChainID = x;
                            break;
                        }
                    }
                    if (iChainID != -1 && drawPos[i] >= 0) {
                        // draw chain name and molecule ID for CGs
                        ig2.drawString(pg.getAllChains().get(iChainID).getPdbChainID(), pl.getFooterStart().x + (drawPos[i] * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));
                        ig2.drawString(pg.getAllChains().get(iChainID).getMacromolID(), pl.getFooterStart().x + (drawPos[i] * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4));
                    }
                    
                    // DEBUG: draw residue number
                    //Boolean drawResNumLabel = true;
                    //if(drawResNumLabel) {
                    //    ig2.drawString(pg.getAllChains().get(iChainID).getPdbChainID() + pg.getVertex(i).shortLabel(), pl.getFooterStart().x + (i * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4));                        
                    //}
                }
            }
            
            // draw "CH" and "ML" strings for complex graphs
            if (iChainID != -1) {
                ig2.drawString("CH", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight * 2 + (stringHeight / 4));
                ig2.drawString("ML", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight * 3 + (stringHeight / 4));
            }
            
            if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
                if (iChainID != -1) {
                    ProteinGraphDrawer.drawLegend(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + lineHeight * 4 + (stringHeight / 4)), pl, pg);                    
                } else {
                    ProteinGraphDrawer.drawLegend(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + lineHeight * 2 + (stringHeight / 4)), pl, pg);
                }
            }
        }
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());
        DrawResult drawRes = new DrawResult(ig2, roi);
        return drawRes;
    }
    
    /**
     * Assigns a color to a chemical property of an AA, useful for drawing residue representations, e.g., in amino acif graphs or matrices.
     * @param chemProp3 the chem prop 3 string, use constants like AminoAcid.CHEMPROP3_AA_STRING_HYDROPHOBIC
     * @return a color, or null if the supplied chemProp3 string is not assigned to any color
     */
    public static Color getChemProp3Color(String chemProp3) {
        Color c = null;
        if(chemProp3.equals(AminoAcid.CHEMPROP3_AA_STRING_HYDROPHOBIC)) {
            c = Color.BLACK;
        }
        if(chemProp3.equals(AminoAcid.CHEMPROP3_AA_STRING_POLAR_CHARGED)) {
            c = Color.RED;
        }
        if(chemProp3.equals(AminoAcid.CHEMPROP3_AA_STRING_POLAR_UNCHARGED)) {
            c = Color.BLUE;
        }
        return c;
    }
    
    /**
     * Assigns a color to a chemical property of an AA, useful for drawing residue representations, e.g., in amino acid graphs or matrices.
     * @param chemProp5 the chem prop 5 string, use constants like AminoAcid.CHEMPROP5_*
     * @return a color, or null if the supplied chemProp5 string is not assigned to any color
     */
    public static Color getChemProp5Color(String chemProp5) {
        Color c = null;
        if(chemProp5.equals(AminoAcid.CHEMPROP5_AA_STRING_HYDROPHOBIC)) {
            c = Color.BLACK;
        }
        if(chemProp5.equals(AminoAcid.CHEMPROP5_AA_STRING_NEGATIVE_CHARGE)) {
            c = Color.ORANGE;
        }
        if(chemProp5.equals(AminoAcid.CHEMPROP5_AA_STRING_POLAR)) {
            c = Color.BLUE;
        }
        if(chemProp5.equals(AminoAcid.CHEMPROP5_AA_STRING_POSITIVE_CHARGE)) {
            c = Color.MAGENTA;
        }
        if(chemProp5.equals(AminoAcid.CHEMPROP5_AA_STRING_SMALL_APOLAR)) {
            c = Color.GRAY;
        }
        return c;
    }
    
    
    /**
     * Draws the drawable graph image of this graph and returns the DrawResult.
     *
     * @return the DrawResult. You can write this to a file or whatever.
     */
    private static DrawResult drawDrawableGraphG2D(IDrawableGraph pg, Map<Integer, String> vertexMarkings) {
        Integer numVerts = pg.getDrawableVertices().size();
        Boolean bw = false;
        
        String graphType = pg.getPropertyString("graphType");
        String pdbid = pg.getPropertyString("pdbid");
        String chainid = pg.getPropertyString("chainid");
        
        PageLayout pl = new PageLayout(numVerts);
        Position2D vertStart = pl.getVertStart();
        SVGGraphics2D ig2;
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        ig2 = new SVGGraphics2D(document);
        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);
        Font font = pl.getStandardFont();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();
        String proteinHeader = "" + graphType + " protein graph of PDB entry " +pdbid + ", chain " + chainid + "";
        Integer stringHeight = fontMetrics.getAscent();
        String sseNumberSeq;
        String sseNumberGraph;
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            ig2.drawString(proteinHeader, pl.headerStart.x, pl.headerStart.y);
        }
        Shape shape, markShape;
        Arc2D.Double arc, markArc;
        ig2.setStroke(new BasicStroke(2));
        Integer edgeType;
        Integer leftVert;
        Integer rightVert;
        Integer leftVertPosX;
        Integer rightVertPosX;
        Integer arcWidth;
        Integer arcHeight;
        Integer arcTopLeftX;
        Integer arcTopLeftY;
        Integer spacerX;
        Integer spacerY;
        Integer iChainID;
        Integer jChainID;
        
        Color markingBorderColor = Color.GRAY;
        Color markingLabelColor = Color.GRAY;
        int vertexMarkingWidth = 4;
        int vertexMarkingLabelOffsetY = 25;
        int edgeMarkingOffsetX = 0; // the number of pixels the gray edge is shifted on the x axis
        
        for (Integer i = 0; i < pg.getDrawableVertices().size(); i++) {
            for (Integer j = i + 1; j < pg.getDrawableVertices().size(); j++) {
                if (pg.containsEdge(i, j)) {
                  
                    edgeType = SpatRel.stringToInt(pg.getSpatRelOfEdge(i, j));
                    if (edgeType.equals(SpatRel.PARALLEL)) {
                        ig2.setPaint(Color.RED);
                    } else if (edgeType.equals(SpatRel.ANTIPARALLEL)) {
                        ig2.setPaint(Color.BLUE);
                    } else if (edgeType.equals(SpatRel.MIXED)) {
                        ig2.setPaint(Color.GREEN);
                    } else if (edgeType.equals(SpatRel.LIGAND)) {
                        ig2.setPaint(Color.MAGENTA);
                    } else if (edgeType.equals(SpatRel.BACKBONE)) {
                        ig2.setPaint(Color.ORANGE);
                    } else if (edgeType.equals(SpatRel.COMPLEX)) {
                        ig2.setPaint(Color.BLACK);
                    } else if (edgeType.equals(SpatRel.OTHER)) {
                        ig2.setPaint(Color.LIGHT_GRAY);
                    } else {
                        ig2.setPaint(Color.LIGHT_GRAY);
                    }
                    if (bw) {
                        ig2.setPaint(Color.LIGHT_GRAY);
                    }
                  
                    if (i < j) {
                        leftVert = i;
                        rightVert = j;
                    } else {
                        leftVert = j;
                        rightVert = i;
                    }
                    leftVertPosX = pl.getVertStart().x + (leftVert * pl.vertDist);
                    rightVertPosX = pl.getVertStart().x + (rightVert * pl.vertDist);
                    arcWidth = rightVertPosX - leftVertPosX;
                    arcHeight = arcWidth / 2;
                    arcTopLeftX = leftVertPosX;
                    arcTopLeftY = pl.getVertStart().y - arcHeight / 2;
                    spacerX = pl.vertRadius;
                    spacerY = 0;
                    
                    Color tmp = ig2.getColor();
                    
                    // mark edges if req. (do this first so real edge gets drawn over this, not vice versa)
                    String vMarking1 = vertexMarkings.get(i);
                    String vMarking2 = vertexMarkings.get(j);
                    if(vMarking1 != null && vMarking2 != null) {
                        ig2.setColor(markingBorderColor);
                        ig2.setStroke(new BasicStroke(6));
                        markArc = new Arc2D.Double(arcTopLeftX + spacerX + edgeMarkingOffsetX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
                        markShape = ig2.getStroke().createStrokedShape(markArc);
                        ig2.fill(markShape);
                    }
                    
                    ig2.setColor(tmp);
                    ig2.setStroke(new BasicStroke(2));
                    
                    arc = new Arc2D.Double(arcTopLeftX + spacerX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
                    shape = ig2.getStroke().createStrokedShape(arc);
                    ig2.fill(shape);                                                                                
                }
            }
        }
        Ellipse2D.Double circle; Ellipse2D.Double markingCircle;
        Rectangle2D.Double rect; Rectangle2D.Double markingRect;
        ig2.setStroke(new BasicStroke(2));
        Integer sseClass;
        for (Integer i = 0; i < pg.getDrawableVertices().size(); i++) {
            sseClass = SSE.sseClassFromFgNotation(pg.getDrawableVertices().get(i).getSseFgNotation());
            
            if (sseClass.equals(SSE.SSECLASS_HELIX)) {
                ig2.setPaint(Color.RED);
            } else if (sseClass.equals(SSE.SSECLASS_BETASTRAND)) {
                ig2.setPaint(Color.BLACK);
            } else if (sseClass.equals(SSE.SSECLASS_LIGAND)) {
                ig2.setPaint(Color.MAGENTA);
            } else if (sseClass.equals(SSE.SSECLASS_OTHER)) {
                ig2.setPaint(Color.GRAY);
            } else {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            if (bw) {
                ig2.setPaint(Color.GRAY);
            }
            
            String vMarking = vertexMarkings.get(i);            
            
            if (sseClass.equals(SSE.SSECLASS_BETASTRAND)) {
                rect = new Rectangle2D.Double(vertStart.x + (i * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());                                
                ig2.fill(rect);
                if(vMarking != null) {
                    //System.out.println("#####Found marking for strand vertex " + i + ".");
                    markingRect = new Rectangle2D.Double(vertStart.x + (i * pl.vertDist), vertStart.y, pl.getVertDiameter() + (0 * vertexMarkingWidth), pl.getVertDiameter() + (0 * vertexMarkingWidth));
                    ig2.setColor(markingBorderColor);
                    ig2.setStroke(new BasicStroke(vertexMarkingWidth));
                    ig2.draw(markingRect);
                    ig2.setColor(markingLabelColor);
                    ig2.drawString(vMarking, new Double(markingRect.getMinX()).intValue() + (pl.vertRadius / 2), new Double(markingRect.getCenterY()).intValue() + vertexMarkingLabelOffsetY);
                }
            } else if (sseClass.equals(SSE.SSECLASS_LIGAND)) {
                int outerStrokeWidth = 4;
                int halfouterStrokeWidth = outerStrokeWidth / 2;
                circle = new Ellipse2D.Double(vertStart.x + (i * pl.vertDist) + halfouterStrokeWidth, vertStart.y + halfouterStrokeWidth, pl.getVertDiameter() - outerStrokeWidth, pl.getVertDiameter() - outerStrokeWidth);
                ig2.setStroke(new BasicStroke(outerStrokeWidth));
                ig2.draw(circle);
                ig2.setStroke(new BasicStroke(2));
                
                // mark vertices if req.
                if(vMarking != null) {
                    //System.out.println("#####Found marking for ligand vertex " + i + ".");
                    markingCircle = new Ellipse2D.Double(vertStart.x + (i * pl.vertDist), vertStart.y, pl.getVertDiameter() + (0 * vertexMarkingWidth), pl.getVertDiameter() + (0 * vertexMarkingWidth));
                    ig2.setColor(markingBorderColor);
                    ig2.setStroke(new BasicStroke(vertexMarkingWidth));
                    ig2.draw(markingCircle);
                    ig2.setColor(markingLabelColor);
                    ig2.drawString(vMarking, new Double(markingCircle.getMinX()).intValue() + (pl.vertRadius / 2), new Double(markingCircle.getCenterY()).intValue() + vertexMarkingLabelOffsetY);
                }
            } else {
                circle = new Ellipse2D.Double(vertStart.x + (i * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(circle);
                if(vMarking != null) {
                    //System.out.println("#####Found marking for other vertex " + i + ".");
                    markingCircle = new Ellipse2D.Double(vertStart.x + (i * pl.vertDist), vertStart.y, pl.getVertDiameter() + (0 * vertexMarkingWidth), pl.getVertDiameter() + (0 * vertexMarkingWidth));
                    ig2.setColor(markingBorderColor);
                    ig2.setStroke(new BasicStroke(vertexMarkingWidth));
                    ig2.draw(markingCircle);
                    ig2.setColor(markingLabelColor);
                    ig2.drawString(vMarking, new Double(markingCircle.getMinX()).intValue() + (pl.vertRadius / 2), new Double(markingCircle.getCenterY()).intValue() + vertexMarkingLabelOffsetY);
                }
            }
            
        
        }
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        if (!bw) {
            if (pg.getDrawableVertices().size() > 0) {
                ig2.drawString("N", vertStart.x - pl.vertDist, vertStart.y + 20);
                ig2.drawString("C", vertStart.x + pg.getDrawableVertices().size() * pl.vertDist, vertStart.y + 20);
            }
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {
            Integer printNth = 1;
            if (pg.getDrawableVertices().size() > 9) {
                printNth = 1;
            }
            if (pg.getDrawableVertices().size() > 99) {
                printNth = 2;
            }
            if (pg.getDrawableVertices().size() > 999) {
                printNth = 3;
            }
            Integer lineHeight = pl.textLineHeight;
            if (pg.getDrawableVertices().size() > 0) {
                ig2.drawString("PG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                //ig2.drawString("SQ", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            iChainID = -1;
            for (Integer i = 0; i < pg.getDrawableVertices().size(); i++) {
                if ((i + 1) % printNth == 0) {
                    sseNumberGraph = "" + (i + 1);
                    //sseNumberSeq = "" + (pg.getVertex(i).getSSESeqChainNum());
                    stringHeight = fontMetrics.getAscent();
                    ig2.drawString(sseNumberGraph, pl.getFooterStart().x + (i * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (stringHeight / 4));
                    //ig2.drawString(sseNumberSeq, pl.getFooterStart().x + (i * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                    
                }
            }
            if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
                ProteinGraphDrawer.drawLegend(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + lineHeight * 2 + (stringHeight / 4)), pl); 
            }
        }
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());
        DrawResult drawRes = new DrawResult(ig2, roi);
        return drawRes;
    }

    /**
     * Creates an encapsulated postscript (EPS) format image file of the folding graph.
     * See http://xmlgraphics.apache.org/commons/postscript.html for documentation on this.
     * @param outputFile the target file
     * @throws IOException In case of an I/O error
     */
    public static void drawProteinGraphToEncPostscriptFile(File outputFile) throws IOException {
        OutputStream out = new FileOutputStream(outputFile);
        out = new BufferedOutputStream(out);
        try {
            EPSDocumentGraphics2D g2d = new EPSDocumentGraphics2D(false);
            g2d.setGraphicContext(new GraphicContext());
            g2d.setupDocument(out, 400, 200);
            g2d.drawRect(0, 0, 400, 200);
            Graphics2D copy = (Graphics2D) g2d.create();
            int c = 12;
            for (int i = 0; i < c; i++) {
                float f = (i + 1) / (float) c;
                Color col = new Color(0.0F, 1 - f, 0.0F);
                copy.setColor(col);
                copy.fillRect(70, 90, 50, 50);
                copy.rotate(-2 * Math.PI / (double) c, 70, 90);
            }
            copy.dispose();
            g2d.rotate(-0.25);
            g2d.setColor(Color.RED);
            g2d.setFont(new Font("sans-serif", Font.PLAIN, 36));
            g2d.drawString("Hello world!", 140, 140);
            g2d.setColor(Color.RED.darker());
            g2d.setFont(new Font("serif", Font.PLAIN, 36));
            g2d.drawString("Hello world!", 140, 180);
            g2d.finish();
        } finally {
            IOUtils.closeQuietly(out);
        }
    }

    /**
     * Draws the legend for the graph at the given position. This legend is not suitable for SEQ folding graphs, because their edges are different.
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position (x, y) where to start drawing
     * @return the x coordinate in the image where the legend ends (which is the left margin + the legend width).
     * This can be used to determine the minimal width of the total image (it has to be at least this value).
     */
    public static Integer drawLegend(SVGGraphics2D ig2, Position2D startPos, PageLayout pl, SSEGraph g) {
        Boolean drawAll = Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend_always_all");
        ig2.setFont(pl.getLegendFont());
        FontMetrics fontMetrics = ig2.getFontMetrics();
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        Integer spacer = 10;
        Integer pixelPosX = startPos.x;
        Integer vertWidth = pl.getVertDiameter();
        Integer vertOffset = pl.getVertDiameter() / 4 * 3;
        String label;
        if (g.numEdges() > 0) {
            label = "[Edges: ";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsContactTypeParallel() || drawAll) {
            label = "parallel";
            ig2.setPaint(Color.RED);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsContactTypeAntiparallel() || drawAll) {
            label = "antiparallel";
            ig2.setPaint(Color.BLUE);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsContactTypeMixed() || drawAll) {
            label = "mixed";
            ig2.setPaint(Color.GREEN);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsContactTypeLigand() || drawAll) {
            label = "ligand";
            ig2.setPaint(Color.MAGENTA);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.getChainEnds().size() > 0 || drawAll) {
            label = "interchain";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.numEdges() > 0) {
            label = "]";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.numVertices() > 0) {
            label = " [Vertices: ";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeHelix() || drawAll) {
            ProteinGraphDrawer.drawSymbolAlphaHelix(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "helix";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeBetaStrand() || drawAll) {
            ProteinGraphDrawer.drawSymbolBetaStrand(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "strand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeLigand() || drawAll) {
            ProteinGraphDrawer.drawSymbolLigand(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "ligand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.containsSSETypeOther() || drawAll) {
            ProteinGraphDrawer.drawSymbolOtherSSE(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "other";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        if (g.numVertices() > 0) {
            label = "]";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        }
        return pixelPosX;
    }
    
    /**
     * Draws the legend for a protein graph, without asking whether it contains certain SSE types or not.
     * @param ig2
     * @param startPos
     * @param pl
     * @return the x pixel position after drawing (with in pixels)
     */
    public static Integer drawLegend(SVGGraphics2D ig2, Position2D startPos, PageLayout pl) {
        Boolean drawAll = Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend_always_all");
        ig2.setFont(pl.getLegendFont());
        FontMetrics fontMetrics = ig2.getFontMetrics();
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        Integer spacer = 10;
        Integer pixelPosX = startPos.x;
        Integer vertWidth = pl.getVertDiameter();
        Integer vertOffset = pl.getVertDiameter() / 4 * 3;
        String label;
        
            label = "[Edges: ";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
        
            label = "parallel";
            ig2.setPaint(Color.RED);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
        
            label = "antiparallel";
            ig2.setPaint(Color.BLUE);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
        
            label = "mixed";
            ig2.setPaint(Color.GREEN);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            label = "ligand";
            ig2.setPaint(Color.MAGENTA);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            label = "interchain";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            label = "]";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            label = " [Vertices: ";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            ProteinGraphDrawer.drawSymbolAlphaHelix(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "helix";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            ProteinGraphDrawer.drawSymbolBetaStrand(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "strand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            ProteinGraphDrawer.drawSymbolLigand(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "ligand";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            ProteinGraphDrawer.drawSymbolOtherSSE(ig2, new Position2D(pixelPosX, startPos.y - vertOffset), pl);
            pixelPosX += vertWidth + spacer;
            label = "other";
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
            label = "]";
            ig2.setPaint(Color.BLACK);
            ig2.drawString(label, pixelPosX, startPos.y);
            pixelPosX += fontMetrics.stringWidth(label) + spacer;
        
        return pixelPosX;
    }

    /**
     * Draws an ADJ folding graph in all image formats, returns a list of written files.
     * @param baseFilePathNoExt the base file path where to put the image (without dot and file extension)
     * @param drawBlackAndWhite whether to omit colors, only useful for non-protein graphs
     * @param formats an array of type DrawTools.IMAGEFORMAT. Do not include SVG, this will always be drawn anyways
     * @param pnfr
     * @return a map of formats to the corresponding output files written to disk
     */
    public static HashMap<DrawTools.IMAGEFORMAT, String> drawFoldingGraphADJ(String baseFilePathNoExt, Boolean drawBlackAndWhite, DrawTools.IMAGEFORMAT[] formats, PTGLNotationFoldResult pnfr) {
        DrawResult drawRes = ProteinGraphDrawer.drawFoldingGraphADJG2D(pnfr);
        String svgFilePath = baseFilePathNoExt + ".svg";
        HashMap<DrawTools.IMAGEFORMAT, String> resultFilesByFormat = new HashMap<DrawTools.IMAGEFORMAT, String>();
        try {
            DrawTools.writeG2dToSVGFile(svgFilePath, drawRes);
            resultFilesByFormat.put(DrawTools.IMAGEFORMAT.SVG, svgFilePath);
            resultFilesByFormat.putAll(DrawTools.convertSVGFileToOtherFormats(svgFilePath, baseFilePathNoExt, drawRes, formats));
        } catch (IOException ex) {
            DP.getInstance().e("Could not write ADJ folding graph file : '" + ex.getMessage() + "'.");
        }
        if (! (Settings.getBoolean("PTGLgraphComputation_B_silent") || Settings.getBoolean("PTGLgraphComputation_B_only_essential_output"))) {
            StringBuilder sb = new StringBuilder();
            sb.append("          Output ADJ folding graph files: ");
            for (DrawTools.IMAGEFORMAT format : resultFilesByFormat.keySet()) {
                sb.append("(").append(format.toString()).append(" => ").append(resultFilesByFormat.get(format)).append(") ");
            }
            System.out.println(sb.toString());
        }
        return resultFilesByFormat;
    }

    /**
     * Draws the SEQ folding graph image of this graph and returns the DrawResult.
     * If 'nonProteinGraph' is true, this graph is considered a custom (=non-protein) graph and the color coding for vertices and edges is NOT used.
     * In that case, the graph is drawn black and white and the labels for the N- and C-termini are NOT drawn.
     *
     * @param pnfr a folding graph notation result
     * @return the DrawResult. You can write this to a file or whatever.
     */
    private static DrawResult drawFoldingGraphSEQG2D(PTGLNotationFoldResult pnfr) {
        FoldingGraph fg = pnfr.getFoldingGraph();
        SSEGraph pg = fg.getParent();
        Integer startVertexInParent = fg.getMinimalVertexIndexInParentGraph();
        Integer endVertexInParent = fg.getMaximalVertexIndexInParentGraph();
        Integer shiftBack = startVertexInParent;
        List<Integer> fgVertexPosInParent = fg.getVertexIndexListInParentGraph();
        List<Integer> parentVertexPosInFG = pg.computeParentGraphVertexPositionsInFoldingGraph(fgVertexPosInParent, pg.getSize());
        Integer numVerts = endVertexInParent - startVertexInParent + 1;
        Boolean bw = false;
        PageLayout pl = new PageLayout(numVerts);
        Position2D vertStart = pl.getVertStart();
        Integer lineHeight = pl.textLineHeight;
        SVGGraphics2D ig2;
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        ig2 = new SVGGraphics2D(document);
        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);
        Font font = pl.getStandardFont();
        Font fontBold = pl.getStandardFontBold();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();
        
        String proteinHeader = "SEQ " + pg.getGraphType() + " folding graph " + fg.getFoldingGraphFoldName() + " (# " + fg.getFoldingGraphNumber() + ") of PDB entry " + pg.getPdbid() + ", chain " + pg.getChainid() + "";
        proteinHeader += "\n" + Settings.getImagesCreatorVersionLine();
        if (FileParser.getWhichSseInfo().equals("author")){
            proteinHeader += "\n" + AUTHOR_SSE_WARNING;
        }
        
        String notation = "SEQ notation: '" + pnfr.seqNotation + "'";
        Integer stringHeight = fontMetrics.getAscent();
        String sseNumberSeq;
        String sseNumberFoldingGraph;
        String sseNumberProteinGraph;
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            DrawTools.drawStringLineBreaks(ig2, proteinHeader, pl.headerStart.x, pl.headerStart.y);
            if (Settings.getBoolean("PTGLgraphComputation_B_print_notations_on_fg_images")) {
                ig2.drawString(notation, pl.headerStart.x, pl.headerStart.y + lineHeight);
            }
        }
        Integer k;
        Integer l;
        Shape shape;
        Arc2D.Double arc;
        ig2.setStroke(new BasicStroke(2));
        Integer edgeType;
        Integer leftVert;
        Integer rightVert;
        Integer leftVertPosX;
        Integer rightVertPosX;
        Integer arcWidth;
        Integer arcHeight;
        Integer arcTopLeftX;
        Integer arcTopLeftY;
        Integer spacerX;
        Integer spacerY;
        List<Integer> fgVertexIndices = fg.getVertexIndexList();
        Integer j;
        for (Integer i = 0; i < fgVertexIndices.size() - 1; i++) {
            j = fgVertexIndices.get(i + 1);
            k = fgVertexPosInParent.get(i);
            l = fgVertexPosInParent.get(j);
            ig2.setPaint(Color.BLACK);
            if (k < l) {
                leftVert = k;
                rightVert = l;
            } else {
                leftVert = l;
                rightVert = k;
            }
            leftVertPosX = pl.getVertStart().x + ((leftVert - shiftBack) * pl.vertDist);
            rightVertPosX = pl.getVertStart().x + ((rightVert - shiftBack) * pl.vertDist);
            arcWidth = rightVertPosX - leftVertPosX;
            arcHeight = arcWidth / 2;
            arcTopLeftX = leftVertPosX;
            arcTopLeftY = pl.getVertStart().y - arcHeight / 2;
            spacerX = pl.vertRadius;
            spacerY = 0;
            arc = new Arc2D.Double(arcTopLeftX + spacerX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
            shape = ig2.getStroke().createStrokedShape(arc);
            ig2.fill(shape);
        }
        Ellipse2D.Double circle;
        Rectangle2D.Double rect;
        ig2.setStroke(new BasicStroke(2));
        for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
            if (pg.getVertex(i).isHelix()) {
                ig2.setPaint(Color.RED);
            } else if (pg.getVertex(i).isBetaStrand()) {
                ig2.setPaint(Color.BLACK);
            } else if (pg.getVertex(i).isLigandSSE()) {
                ig2.setPaint(Color.MAGENTA);
            } else if (pg.getVertex(i).isOtherSSE()) {
                ig2.setPaint(Color.GRAY);
            } else {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            Integer parentVertexPosInFoldingGraph = parentVertexPosInFG.get(i);
            if (parentVertexPosInFoldingGraph < 0) {
                ig2.setPaint(Color.LIGHT_GRAY);
            }
            if (bw) {
                ig2.setPaint(Color.GRAY);
            }
            SSE sse = pg.getVertex(i);
            if (sse.isBetaStrand()) {
                rect = new Rectangle2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(rect);
            } else if (sse.isLigandSSE()) {
                circle = new Ellipse2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.setStroke(new BasicStroke(3));
                ig2.draw(circle);
                ig2.setStroke(new BasicStroke(2));
            } else {
                circle = new Ellipse2D.Double(vertStart.x + ((i - shiftBack) * pl.vertDist), vertStart.y, pl.getVertDiameter(), pl.getVertDiameter());
                ig2.fill(circle);
            }
        }
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        if (!bw) {
            if (fg.getSize() > 0) {
                ig2.setFont(fontBold);
                ig2.drawString("N", vertStart.x - pl.vertDist, vertStart.y + 20);
                ig2.drawString("C", vertStart.x + numVerts * pl.vertDist, vertStart.y + 20);
                ig2.setFont(font);
            }
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {
            Integer printNth = 1;
            if (fg.getSize() > 99) {
                printNth = 2;
            }
            if (fg.getSize() > 999) {
                printNth = 3;
            }
            if (fg.getSize() > 0) {
                ig2.drawString("FG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                ig2.setPaint(Color.LIGHT_GRAY);
                ig2.drawString("SQ", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                ig2.drawString("PG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + lineHeight + (stringHeight / 4));
                ig2.setPaint(Color.BLACK);
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            ig2.setPaint(Color.BLACK);
            for (Integer i = startVertexInParent; i <= endVertexInParent; i++) {
                if ((i + 1) % printNth == 0) {
                    sseNumberFoldingGraph = "" + (parentVertexPosInFG.get(i) >= 0 ? (parentVertexPosInFG.get(i) + 1) : "");
                    sseNumberSeq = "" + (pg.getVertex(i).getSSESeqChainNum());
                    sseNumberProteinGraph = "" + (i + 1);
                    stringHeight = fontMetrics.getAscent();
                    ig2.drawString(sseNumberFoldingGraph, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (stringHeight / 4));
                    ig2.setPaint(Color.LIGHT_GRAY);
                    ig2.drawString(sseNumberSeq, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                    ig2.drawString(sseNumberProteinGraph, pl.getFooterStart().x + ((i - shiftBack) * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));
                    ig2.setPaint(Color.BLACK);
                }
            }
            if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
                ProteinGraphDrawer.drawLegendSEQ(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4)), pl, fg);
            }
        }
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());
        DrawResult drawRes = new DrawResult(ig2, roi);
        return drawRes;
    }

    /**
     * Draws a KEY folding graph in all formats, returns a list of written files.
     * @param baseFilePathNoExt the base file path where to put the image (without dot and file extension)
     * @param drawBlackAndWhite whether to omit colors, only useful for non-protein graphs
     * @param formats an array of type DrawTools.IMAGEFORMAT. Do not include SVG, this will always be drawn anyways
     * @param pnfr
     * @return a map of formats to the corresponding output files written to disk
     */
    public static HashMap<DrawTools.IMAGEFORMAT, String> drawFoldingGraphKEY(String baseFilePathNoExt, Boolean drawBlackAndWhite, DrawTools.IMAGEFORMAT[] formats, PTGLNotationFoldResult pnfr) {
        HashMap<DrawTools.IMAGEFORMAT, String> resultFilesByFormat = new HashMap<DrawTools.IMAGEFORMAT, String>();
        if (!pnfr.getFoldingGraph().supportsKeyNotation()) {
            if (pnfr.keyNotation != null || pnfr.keyStartFG >= 0) {
                DP.getInstance().e("SSEGraph", "drawFoldingGraphKEY: FG reported to NOT support KEY notation, but a KEY notation string is set: '" + pnfr.keyNotation + "', start index is parent is " + pnfr.keyStartFG + ".");
            }
            if (! (Settings.getBoolean("PTGLgraphComputation_B_silent") || Settings.getBoolean("PTGLgraphComputation_B_only_essential_output"))) {
                System.out.println("          Output KEY folding graph files: none, FG does not support KEY notation.");
            }
            return resultFilesByFormat;
        }
        DrawResult drawRes = ProteinGraphDrawer.drawFoldingGraphKEYG2D(pnfr);
        if (drawRes == null) {
            return resultFilesByFormat;
        }
        String svgFilePath = baseFilePathNoExt + ".svg";
        try {
            DrawTools.writeG2dToSVGFile(svgFilePath, drawRes);
            resultFilesByFormat.put(DrawTools.IMAGEFORMAT.SVG, svgFilePath);
            resultFilesByFormat.putAll(DrawTools.convertSVGFileToOtherFormats(svgFilePath, baseFilePathNoExt, drawRes, formats));
        } catch (IOException ex) {
            DP.getInstance().e("Could not write KEY folding graph file : '" + ex.getMessage() + "'.");
        }
        if (! (Settings.getBoolean("PTGLgraphComputation_B_silent") || Settings.getBoolean("PTGLgraphComputation_B_only_essential_output"))) {
            StringBuilder sb = new StringBuilder();
            sb.append("          Output KEY folding graph files: ");
            for (DrawTools.IMAGEFORMAT format : resultFilesByFormat.keySet()) {
                sb.append("(").append(format.toString()).append(" => ").append(resultFilesByFormat.get(format)).append(") ");
            }
            System.out.println(sb.toString());
        }
        return resultFilesByFormat;
    }

    /**
     * Draws the KEY folding graph image of this graph and returns the DrawResult. This is the new version, started JAN 2015.
     * If 'nonProteinGraph' is true, this graph is considered a custom (=non-protein) graph and the color coding for vertices and edges is NOT used.
     * In that case, the graph is drawn black and white and the labels for the N- and C-termini are NOT drawn.
     *
     * See http://ptgl.uni-frankfurt.de/cgi-bin/showpict.pl?topology=a&rep=3&protlist=7timA+7timB+&nmrlist=
     * and http://ptgl.uni-frankfurt.de/ptglhelp.html#key
     *
     * See the KEY beta FG #3 of 1GOS chain A as an example.
     *
     * @param pnfr a folding graph notation result
     * @return the DrawResult. You can write this to a file or whatever.
     */
    private static DrawResult drawFoldingGraphKEYG2D(PTGLNotationFoldResult pnfr) {
        FoldingGraph fg = pnfr.getFoldingGraph();
        SSEGraph pg = fg.getParent();
        boolean debug = false;
        /*
        if(pg.pdbid.equals("8icd") && pg.chainid.equals("A")) {
        if(fg.graphType.equals("beta") && fg.getFoldingGraphNumber().equals(1)) {
        debug = true;
        }
        }
         */
        if (debug) {
            System.out.println("******************** START *********************");
        }
        if (fg.isBifurcated()) {
            DP.getInstance().e("SSEGraph", "drawFoldingGraphKEYG2D: This FG is bifurcated, KEY notation not supported.");
            return null;
        }
        Integer leftMostVertexInParent = fg.getMinimalVertexIndexInParentGraph();
        Integer rightMostVertexInParent = fg.getMaximalVertexIndexInParentGraph();
        List<Integer> keyposParentIndicesSeqOrder = pnfr.keypos;
        Integer keystartFGIndex = pnfr.keyStartFG;
        if (keystartFGIndex == null) {
            DP.getInstance().e("SSEGraph", "drawFoldingGraphKEYG2D: keystartFGIndex is null");
        }
        fg.computeSpatialVertexOrdering();
        List<Integer> keyposFGIndicesSpatOrder = fg.getSpatOrder();
        if (keyposFGIndicesSpatOrder.size() != fg.getSize()) {
            if (fg.isASingleCycle()) {
                keyposFGIndicesSpatOrder = fg.getSpatialOrderingOfVertexIndicesForSingleCycleFG(keystartFGIndex);
                fg.setSpatOrder(keyposFGIndicesSpatOrder);
                if (keyposFGIndicesSpatOrder.size() != fg.getSize()) {
                    DP.getInstance().e("SSEGraph", "Could not draw KEY notation: Spatorder size " + keyposFGIndicesSpatOrder.size() + " of chain " + fg.getChainid() + " gt " + fg.getGraphType() + " FG " + fg.getFoldingGraphFoldName() + " (#" + fg.getFoldingGraphNumber() + ") does not match FG size " + fg.getSize() + " even when allowing circles. Parent verts of FG: " + IO.intListToString(fg.getVertexIndexListInParentGraph()) + ". KEY='" + pnfr.keyNotation + "'. " + (fg.isASingleCycle() ? "FG is a single cycle." : "FG is NOT a single cycle."));
                    return null;
                } else {
                    //System.err.println("Single cycle handling fixed it.");
                }
            } else {
                DP.getInstance().e("SSEGraph", "Could not draw KEY notation: Spatorder size " + keyposFGIndicesSpatOrder.size() + " of chain " + fg.getChainid() + " gt " + fg.getGraphType() + " FG " + fg.getFoldingGraphFoldName() + " (#" + fg.getFoldingGraphNumber() + ") does not match FG size " + fg.getSize() + ". Parent verts of FG: " + IO.intListToString(fg.getVertexIndexListInParentGraph()) + ". KEY='" + pnfr.keyNotation + "'. " + (fg.isASingleCycle() ? "FG is a single cycle." : "FG is NOT a single cycle."));
                return null;
            }
        }
        Integer keyendFGIndex = keyposFGIndicesSpatOrder.get(fg.getSize() - 1);
        
        if (debug) {
            DP.getInstance().d("keyposFGIndicesSpatOrder=" + IO.intListToString(keyposFGIndicesSpatOrder));
        }
        List<Integer> spatOrderseqDistToPrevious = new ArrayList<>();
        spatOrderseqDistToPrevious.add(null);
        Integer cur;
        Integer last;
        Integer spatPosCur;
        Integer spatPosLast;
        for (Integer i = 1; i < fg.getSize(); i++) {
            cur = i;
            last = i - 1;
            spatPosCur = keyposFGIndicesSpatOrder.indexOf(cur);
            spatPosLast = keyposFGIndicesSpatOrder.indexOf(last);
            spatOrderseqDistToPrevious.add(spatPosCur - spatPosLast);
        }
        List<Integer> fgVertexPosInParent = fg.getVertexIndexListInParentGraph();
        List<Integer> parentVertexPosInFG = ProtGraph.computeParentGraphVertexPositionsInFoldingGraph(fgVertexPosInParent, pg.getSize());
        Integer numVerts = fg.getSize();
        Boolean bw = false;
        PageLayout pl = new PageLayout(numVerts);
        pl.isForKEY = true;
        Position2D vertStart = pl.getVertStart();
        Integer lineHeight = pl.textLineHeight;
        Integer vertDist = 50;
        Integer vertHeight = 80;
        Integer vertWidth = 40;
        Integer vertStartX = pl.getVertStart().x;
        Integer vertStartY = pl.getVertStart().y;
        SVGGraphics2D ig2;
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        ig2 = new SVGGraphics2D(document);
        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);
        Font font = pl.getStandardFont();
        Font fontBold = pl.getStandardFontBold();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();
        
        String proteinHeader = "KEY " + pg.getGraphType() + " folding graph " + fg.getFoldingGraphFoldName() + " (# " + fg.getFoldingGraphNumber() + ") of PDB entry " + pg.getPdbid() + ", chain " + pg.getChainid() + "";
        proteinHeader += "\n" + Settings.getImagesCreatorVersionLine();
        if (FileParser.getWhichSseInfo().equals("author")){
            proteinHeader += "\n" + AUTHOR_SSE_WARNING;
        }
        
        String notation = "KEY notation: '" + pnfr.keyNotation + "'";
        if (debug) {
            DP.getInstance().d("Notation from PNFR: " + notation);
        }
        Integer stringHeight = fontMetrics.getAscent();
        String sseNumberSeqInChain;
        String sseNumberFoldingGraph;
        String sseNumberProteinGraph;
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            DrawTools.drawStringLineBreaks(ig2, proteinHeader, pl.headerStart.x, pl.headerStart.y);
        }
        Integer[] shiftBack = new Integer[pg.getSize()];
        Integer shift = 0;
        for (int i = 0; i < parentVertexPosInFG.size(); i++) {
            if (parentVertexPosInFG.get(i) < 0) {
                shift++;
            }
            shiftBack[i] = shift;
        }
        boolean debugDrawingKEY = false;
        //if(fg.toShortString().equals("1gos-B-alpha-FG0[4V,3E]")) {
        Integer[] orientationsSpatOrder = new Integer[fg.getSize()];
        Integer[] orientationsSeqOrder = new Integer[fg.getSize()];
        List<Integer> keyposFGIndices = new ArrayList<>();
        Integer v;
        for (int i = 0; i < keyposParentIndicesSeqOrder.size(); i++) {
            v = parentVertexPosInFG.get(keyposParentIndicesSeqOrder.get(i));
            if (v < 0) {
                DP.getInstance().e("SSEGRaph", "Keypos parent index #" + i + " is " + v + ", skipping.");
            } else {
                keyposFGIndices.add(v);
            }
        }
        if (keyposFGIndices.size() != fg.getSize()) {
            DP.getInstance().e("SSEGRaph", "keyposFGIndices.size()=" + keyposFGIndices.size() + ", fg.size=" + fg.getSize() + ".");
        }
        Integer firstVertexSpatFGIndex = fg.getSpatOrder().get(0);
        orientationsSpatOrder[0] = FoldingGraph.ORIENTATION_UPWARDS;
        orientationsSeqOrder[firstVertexSpatFGIndex] = FoldingGraph.ORIENTATION_UPWARDS;
        if (debugDrawingKEY) {
            System.out.println("Setting orientationsSeqOrder[" + firstVertexSpatFGIndex + "] to " + FoldingGraph.ORIENTATION_UPWARDS + ".");
        }
        Integer currentVert;
        Integer lastVert;
        StringBuilder KEYNotation = new StringBuilder();
        Boolean appendSSEType = true;
        String fgGraphtType = pnfr.getFoldingGraph().getGraphType();
        if (fgGraphtType.equals(SSEGraph.GRAPHTYPE_ALPHA) || fgGraphtType.equals(SSEGraph.GRAPHTYPE_BETA)) {
        }
        String bracketStart = "[";
        String bracketEnd = "]";
        if (fg.isBifurcated()) {
            bracketStart = "{";
            bracketEnd = "}";
        } else if (fg.isASingleCycle()) {
            bracketStart = "(";
            bracketEnd = ")";
        }
        KEYNotation.append(bracketStart);
        if (appendSSEType) {
            KEYNotation.append(fg.getVertex(keystartFGIndex).getLinearNotationLabel());
        }
        if (keyposFGIndices.size() > 1) {
            for (int i = 1; i < keyposFGIndicesSpatOrder.size(); i++) {
                if (i < (keyposFGIndicesSpatOrder.size())) {
                    if (!appendSSEType && i == 1) {
                    } else {
                        KEYNotation.append(",");
                    }
                }
                currentVert = keyposFGIndicesSpatOrder.get(i);
                lastVert = keyposFGIndicesSpatOrder.get(i - 1);
                KEYNotation.append(currentVert - lastVert);
                Integer spatRel = fg.getContactSpatRel(lastVert, currentVert);
                if (Objects.equals(spatRel, SpatRel.PARALLEL)) {
                    KEYNotation.append("x");
                    orientationsSpatOrder[i] = orientationsSpatOrder[i - 1];
                    orientationsSeqOrder[currentVert] = orientationsSpatOrder[i];
                } else {
                    orientationsSpatOrder[i] = (Objects.equals(orientationsSpatOrder[i - 1], FoldingGraph.ORIENTATION_UPWARDS) ? FoldingGraph.ORIENTATION_DOWNWARDS : FoldingGraph.ORIENTATION_UPWARDS);
                    orientationsSeqOrder[currentVert] = orientationsSpatOrder[i];
                }
                if (appendSSEType) {
                    KEYNotation.append(fg.getVertex(currentVert).getLinearNotationLabel());
                }
            }
        }
        KEYNotation.append(bracketEnd);
        if (debug) {
            DP.getInstance().d("Notation from this func: " + KEYNotation.toString());
        }
        Integer[] newOrientations = new Integer[fg.getSpatOrder().size()];
        Arrays.fill(newOrientations, FoldingGraph.ORIENTATION_NONE);
        Integer currentVertexIndexInFGSequential;
        List<Shape> connShapes;
        Polygon pol;
        Position2D p;
        Position2D lastP;
        for (int i = 0; i < fg.getSpatOrder().size(); i++) {
            currentVertexIndexInFGSequential = fg.getSpatOrder().get(i);
            p = new Position2D(vertStartX + (i * vertDist) + pl.vertRadius / 2, vertStartY);
            Integer contactTypeInt = null;
            if (i > 0) {
                contactTypeInt = fg.getContactSpatRel(fg.getSpatOrder().get(i - 1), fg.getSpatOrder().get(i));
                if (contactTypeInt.equals(SpatRel.PARALLEL)) {
                    newOrientations[i] = newOrientations[i - 1];
                } else {
                    newOrientations[i] = (Objects.equals(newOrientations[i - 1], FoldingGraph.ORIENTATION_UPWARDS) ? FoldingGraph.ORIENTATION_DOWNWARDS : FoldingGraph.ORIENTATION_UPWARDS);
                }
            } else {
                newOrientations[i] = FoldingGraph.ORIENTATION_UPWARDS;
            }
            ig2.setColor(Color.BLACK);
            if (fg.getVertex(currentVertexIndexInFGSequential).isBetaStrand()) {
                ig2.setColor(Color.BLACK);
                pol = DrawTools.getDefaultArrowPolygonLowestPointAt(p.x, p.y, newOrientations[i]);
            } else if (fg.getVertex(currentVertexIndexInFGSequential).isHelix()) {
                ig2.setColor(Color.RED);
                pol = DrawTools.getDefaultBarrelPolygonLowestPointAt(p.x, p.y);
            } else if (fg.getVertex(currentVertexIndexInFGSequential).isLigandSSE()) {
                ig2.setColor(Color.MAGENTA);
                pol = DrawTools.getDefaultBarrelPolygonLowestPointAt(p.x, p.y);
            } else {
                ig2.setColor(Color.LIGHT_GRAY);
                pol = DrawTools.getDefaultBarrelPolygonLowestPointAt(p.x, p.y);
            }
            ig2.draw(pol);
            ig2.fill(pol);
            ig2.setColor(Color.BLACK);
            ig2.drawString((currentVertexIndexInFGSequential + 1) + "", pl.getFooterStart().x + (i * vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (stringHeight / 4));
            List<Integer> parentMapping = fg.getVertexIndexListInParentGraph();
            sseNumberProteinGraph = "" + (parentMapping.get(currentVertexIndexInFGSequential) + 1);
            sseNumberSeqInChain = "" + (pg.getSSEBySeqPosition(parentMapping.get(currentVertexIndexInFGSequential)).getSSESeqChainNum());
            ig2.setColor(Color.LIGHT_GRAY);
            ig2.drawString(sseNumberSeqInChain, pl.getFooterStart().x + (i * vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
            ig2.drawString(sseNumberProteinGraph, pl.getFooterStart().x + (i * vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));
            ig2.setColor(Color.BLACK);
        }
        if (debug) {
            System.out.print("Orientations: ");
            for (Integer o : newOrientations) {
                System.out.print(FoldingGraph.getOrientationString(o) + " ");
            }
            System.out.println("");
        }
        StringBuilder testKEY = new StringBuilder();
        String[] newKey = new String[fg.getSize() - 1];
        Arrays.fill(newKey, "?");
        ig2.setColor(Color.BLACK);
        for (int i = 0; i < fg.getSpatOrder().size(); i++) {
            currentVertexIndexInFGSequential = fg.getSpatOrder().get(i);
            p = new Position2D(vertStartX + (i * vertDist) + pl.vertRadius / 2, vertStartY);
            Integer previousVertexIndexInFGSequential;
            Integer previousVertexIndexSpatial;
            Integer currentVertexIndexSpatial;
            Integer relDist;
            String relDistString;
            List<Shape> shapes;
            if (currentVertexIndexInFGSequential > 0) {
                previousVertexIndexInFGSequential = currentVertexIndexInFGSequential - 1;
                currentVertexIndexSpatial = fg.getSpatOrder().indexOf(currentVertexIndexInFGSequential);
                previousVertexIndexSpatial = fg.getSpatOrder().indexOf(previousVertexIndexInFGSequential);
                relDist = currentVertexIndexSpatial - previousVertexIndexSpatial;
                relDistString = relDist + "";
                if (appendSSEType) {
                    relDistString += fg.getVertex(currentVertexIndexInFGSequential).getLinearNotationLabel();
                }
                newKey[previousVertexIndexInFGSequential] = relDistString;
                if (debug) {
                    System.out.println("EDGE (i=" + i + ") This is vert " + currentVertexIndexInFGSequential + " at position " + currentVertexIndexSpatial + ", last vert sequential was " + previousVertexIndexInFGSequential + " at position " + previousVertexIndexSpatial + ".");
                }
                lastP = new Position2D(vertStartX + (previousVertexIndexSpatial * vertDist) + pl.vertRadius / 2, vertStartY);
                Integer lastOrientation = newOrientations[previousVertexIndexSpatial];
                Integer currentOrientation = newOrientations[currentVertexIndexSpatial];
                if (debug) {
                    System.out.println("  lastOrientation was: " + FoldingGraph.getOrientationString(lastOrientation) + ", current is " + FoldingGraph.getOrientationString(currentOrientation) + ". relDist(" + previousVertexIndexSpatial + " to " + currentVertexIndexSpatial + ")=" + relDist + ".");
                }
                if (lastOrientation.equals(currentOrientation)) {
                    if (lastOrientation.equals(FoldingGraph.ORIENTATION_UPWARDS)) {
                        if (debug) {
                            ig2.setColor(Color.RED);
                        }
                        if (Settings.getBoolean("PTGLgraphComputation_B_key_use_alternate_arcs")) {
                            shapes = DrawTools.getCrossoverArcConnectorAlternativeBezierVersion(lastP.x, lastP.y - vertHeight, p.x, p.y, ig2.getStroke(), true, 0);
                        } else {
                            shapes = DrawTools.getCrossoverArcConnector(lastP.x, lastP.y - vertHeight, p.x, p.y, ig2.getStroke(), true, 0);
                        }
                    } else {
                        if (debug) {
                            ig2.setColor(Color.ORANGE);
                        }
                        if (Settings.getBoolean("PTGLgraphComputation_B_key_use_alternate_arcs")) {
                            shapes = DrawTools.getCrossoverArcConnectorAlternativeBezierVersion(lastP.x, lastP.y, p.x, p.y - vertHeight, ig2.getStroke(), false, 0);
                        } else {
                            shapes = DrawTools.getCrossoverArcConnector(lastP.x, lastP.y, p.x, p.y - vertHeight, ig2.getStroke(), false, 0);
                        }
                    }
                } else {
                    if (currentOrientation.equals(FoldingGraph.ORIENTATION_UPWARDS)) {
                        if (debug) {
                            ig2.setColor(Color.BLUE);
                        }
                        shapes = DrawTools.getArcConnector(lastP.x, lastP.y, p.x, p.y, ig2.getStroke(), false, 0);
                    } else {
                        if (debug) {
                            ig2.setColor(Color.GREEN);
                        }
                        shapes = DrawTools.getArcConnector(lastP.x, lastP.y - vertHeight, p.x, p.y - vertHeight, ig2.getStroke(), true, 0);
                    }
                }
                for (Shape s : shapes) {
                    ig2.draw(s);
                }
                ig2.setColor(Color.BLACK);
            }
            Integer terminusCurrentVertexIndexSpatial = fg.getSpatOrder().indexOf(currentVertexIndexInFGSequential);
            Integer terminusCurrentOrientation = newOrientations[terminusCurrentVertexIndexSpatial];
            ig2.setFont(fontBold);
            if (currentVertexIndexInFGSequential.equals(0)) {
                ig2.setColor(Color.BLACK);
                ig2.drawString("N", p.x - 5, p.y + 25);
            }
            if (currentVertexIndexInFGSequential.equals(fg.getSpatOrder().size() - 1)) {
                ig2.setColor(Color.BLACK);
                if (terminusCurrentOrientation.equals(FoldingGraph.ORIENTATION_UPWARDS)) {
                    ig2.drawString("C", p.x - 5, p.y - 90);
                } else {
                    ig2.drawString("C", p.x - 5, p.y + 25);
                }
            }
            ig2.setFont(font);
        }
        testKEY.append(bracketStart);
        if (appendSSEType) {
            testKEY.append(fg.getVertex(fg.getSpatOrder().get(0)).getLinearNotationLabel());
            testKEY.append(",");
        }
        for (int i = 0; i < newKey.length; i++) {
            testKEY.append(newKey[i]);
            if (i < (newKey.length - 1)) {
                testKEY.append(",");
            }
        }
        testKEY.append(bracketEnd);
        if (debug) {
            System.out.println("newKEY = " + testKEY.toString());
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            if (Settings.getBoolean("PTGLgraphComputation_B_print_notations_on_fg_images")) {
                ig2.drawString(notation, pl.headerStart.x, pl.headerStart.y + lineHeight);
            }
        }
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {
            if (fg.getSize() > 0) {
                ig2.setPaint(Color.BLACK);
                ig2.drawString("FG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                ig2.setPaint(Color.LIGHT_GRAY);
                ig2.drawString("SQ", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                ig2.drawString("PG", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + lineHeight + (stringHeight / 4));
                ig2.setPaint(Color.BLACK);
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
                ProteinGraphDrawer.drawLegendKEY(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4)), pl, fg);
            }
        }
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());
        DrawResult drawRes = new DrawResult(ig2, roi);
        if (debug) {
            System.out.println("******************** END *********************");
        }
        return drawRes;
    }

    /**
     * Draws an DEF folding graph in all image formats, returns a list of written files.
     * @param baseFilePathNoExt the base file path where to put the image (without dot and file extension)
     * @param drawBlackAndWhite whether to omit colors, only useful for non-protein graphs
     * @param formats an array of type DrawTools.IMAGEFORMAT. Do not include SVG, this will always be drawn anyways
     * @param pnfr
     * @return a map of formats to the corresponding output files written to disk
     */
    public static HashMap<DrawTools.IMAGEFORMAT, String> drawFoldingGraphDEF(String baseFilePathNoExt, Boolean drawBlackAndWhite, DrawTools.IMAGEFORMAT[] formats, PTGLNotationFoldResult pnfr) {
        DrawResult drawRes = ProteinGraphDrawer.drawFoldingGraphDEFG2D(pnfr);
        String svgFilePath = baseFilePathNoExt + ".svg";
        HashMap<DrawTools.IMAGEFORMAT, String> resultFilesByFormat = new HashMap<DrawTools.IMAGEFORMAT, String>();
        try {
            DrawTools.writeG2dToSVGFile(svgFilePath, drawRes);
            resultFilesByFormat.put(DrawTools.IMAGEFORMAT.SVG, svgFilePath);
            resultFilesByFormat.putAll(DrawTools.convertSVGFileToOtherFormats(svgFilePath, baseFilePathNoExt, drawRes, formats));
        } catch (IOException ex) {
            DP.getInstance().e("Could not write DEF folding graph file : '" + ex.getMessage() + "'.");
        }
        if (! (Settings.getBoolean("PTGLgraphComputation_B_silent") || Settings.getBoolean("PTGLgraphComputation_B_only_essential_output"))) {
            StringBuilder sb = new StringBuilder();
            sb.append("          Output DEF folding graph files: ");
            for (DrawTools.IMAGEFORMAT format : resultFilesByFormat.keySet()) {
                sb.append("(").append(format.toString()).append(" => ").append(resultFilesByFormat.get(format)).append(") ");
            }
            System.out.println(sb.toString());
        }
        return resultFilesByFormat;
    }

    /**
     * Draws the symbol for SSEs of type other at the given position.
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position where to draw
     * @param pl the PageLayout to use (determines the width and height)
     */
    protected static void drawSymbolOtherSSE(SVGGraphics2D ig2, Position2D startPos, PageLayout pl) {
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.GRAY);
        Ellipse2D.Double circle = new Ellipse2D.Double(startPos.x, startPos.y, pl.getVertDiameter(), pl.getVertDiameter());
        ig2.fill(circle);
    }

    /**
     * Draws the 3 symbols for parent graph SSEs (gray alpha/strand/ligand).
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position where to draw
     * @param pl the PageLayout to use (determines the width and height)
     */
    protected static void drawSymbolsParentSSEs(SVGGraphics2D ig2, Position2D startPos, PageLayout pl) {
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.GRAY);
        Ellipse2D.Double circle = new Ellipse2D.Double(startPos.x, startPos.y, pl.getVertDiameter(), pl.getVertDiameter());
        ig2.fill(circle);
        Rectangle2D.Double rect = new Rectangle2D.Double(startPos.x + 25, startPos.y, pl.getVertDiameter(), pl.getVertDiameter());
        ig2.fill(rect);
        ig2.setStroke(new BasicStroke(3));
        Ellipse2D.Double circle2 = new Ellipse2D.Double(startPos.x + 50, startPos.y, pl.getVertDiameter(), pl.getVertDiameter());
        ig2.draw(circle2);
        ig2.setStroke(new BasicStroke(2));
    }

    /**
     * Draws the symbol for a ligand at the given position.
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position where to draw
     * @param pl the PageLayout to use (determines the width and height)
     */
    protected static void drawSymbolLigand(SVGGraphics2D ig2, Position2D startPos, PageLayout pl) {
        ig2.setStroke(new BasicStroke(3));
        ig2.setPaint(Color.MAGENTA);
        Ellipse2D.Double circle = new Ellipse2D.Double(startPos.x, startPos.y, pl.getVertDiameter(), pl.getVertDiameter());
        ig2.draw(circle);
        ig2.setStroke(new BasicStroke(2));
    }

    /**
     * Draws the symbol for an alpha helix at the given position.
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position where to draw
     * @param pl the PageLayout to use (determines the width and height)
     */
    protected static void drawSymbolAlphaHelix(SVGGraphics2D ig2, Position2D startPos, PageLayout pl) {
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.RED);
        Ellipse2D.Double circle = new Ellipse2D.Double(startPos.x, startPos.y, pl.getVertDiameter(), pl.getVertDiameter());
        ig2.fill(circle);
    }

    /**
     * Draws the symbol for a beta strand at the given position.
     * @param ig2 the SVGGraphics2D object on which to draw
     * @param startPos the start position where to draw
     * @param pl the PageLayout to use (determines the width and height)
     */
    protected static void drawSymbolBetaStrand(SVGGraphics2D ig2, Position2D startPos, PageLayout pl) {
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);
        Rectangle2D.Double rect = new Rectangle2D.Double(startPos.x, startPos.y, pl.getVertDiameter(), pl.getVertDiameter());
        ig2.fill(rect);
    }
    
}
