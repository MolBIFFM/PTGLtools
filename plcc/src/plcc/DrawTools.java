/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package plcc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.batik.apps.rasterizer.DestinationType;
import org.apache.batik.apps.rasterizer.SVGConverter;
import org.apache.batik.apps.rasterizer.SVGConverterException;
import org.apache.batik.svggen.SVGGraphics2DIOException;
import tools.DP;

/**
 *
 * @author ts
 */
public class DrawTools {
    
    public static final String FORMAT_PNG = "PNG";
    public static final String FORMAT_SVG = "SVG";
    public static final String FORMAT_TIFF = "TIFF";
    public static final String FORMAT_PDF = "PDF";
    public static final String FORMAT_JPEG = "JPEG";
    
    public enum IMAGEFORMAT { PNG, SVG, TIFF, PDF, JPEG }
    
    
    /**
     * Writes an SVGGraphics object (the contents of the canvas) to an image file in SVG format.
     * @param svgFilePath the path to the output file, including file extension
     * @param drawRes the graphics object, its contents will be saved to the image file
     * @throws IOException if something went wrong with writing the file
     */
    public static void writeG2dToSVGFile(String svgFilePath, DrawResult drawRes) throws IOException {
        drawRes.g2d.stream(new FileWriter(svgFilePath), false);                                                    
    }
    
    /**
     * Converts the input SVG file to various other formats.
     * @param svgInputFilePath the SVG input image
     * @param outputFileBasePathNoExt the base output file name (without . and without extension)
     * @param drawRes the drawRes, which is required to determine the ROI within the SVG
     * @param formats a list of formats, use the constants in DrawTools class
     * @return a list of files that were written successfully, by format
     */
    public static HashMap<IMAGEFORMAT, String> convertSVGFileToOtherFormats(String svgInputFilePath, String outputFileBasePathNoExt, DrawResult drawRes, IMAGEFORMAT[] formats) {
        
        HashMap<IMAGEFORMAT, String> outfilesByFormat = new HashMap<IMAGEFORMAT, String>();
        
        // write other formats
        SVGConverter svgConverter = new SVGConverter();
        svgConverter.setArea(drawRes.roi);
        svgConverter.setWidth((float) drawRes.roi.getWidth());
        svgConverter.setHeight((float) drawRes.roi.getHeight());

        String formatFileExt = "";
        for(IMAGEFORMAT format : formats) {           

            if(format.equals(IMAGEFORMAT.PNG)) {                
                svgConverter.setDestinationType(DestinationType.PNG);                
                formatFileExt = DestinationType.PNG_EXTENSION;
            } else if(format.equals(IMAGEFORMAT.JPEG)) {
                svgConverter.setDestinationType(DestinationType.JPEG);
                svgConverter.setQuality(0.8F);  // JPEG compression
                formatFileExt = DestinationType.JPEG_EXTENSION;
            } else if(format.equals(IMAGEFORMAT.TIFF)) {
                svgConverter.setDestinationType(DestinationType.TIFF);                
                formatFileExt = DestinationType.TIFF_EXTENSION;
            } else if(format.equals(IMAGEFORMAT.PDF)) {
                svgConverter.setDestinationType(DestinationType.PDF);
                formatFileExt = DestinationType.PDF_EXTENSION;
            } else {
                DP.getInstance().w("Unsupported image output format ignored.");
                continue;
            }

            svgConverter.setSources(new String[]{svgInputFilePath});
            svgConverter.setDst(new File(outputFileBasePathNoExt));
            try {      
                svgConverter.execute();
                outfilesByFormat.put(format, outputFileBasePathNoExt + formatFileExt);
            } catch (SVGConverterException ex) {
                DP.getInstance().e("Could not convert SVG file to format '" + format + "', skipping.");
            }
        }
        return outfilesByFormat;
    }
    
}