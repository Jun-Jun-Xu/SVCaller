/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import app.*;
import app.Deletion;
import app.Insertion.*;
import app.genotypeDel.*;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import pgl.infra.utils.CLIInterface;


public class Entrance implements CLIInterface {
    Options options = new Options();
    String introduction = this.createIntroduction();
    String app = null;

    String inputFile = null;
    String outputFile=null;
    String length = null ;
    String quality = null ;
    String path = null;
    String parameterFiles = null;
    String threads = null;
    String chrNum = null;
    String range = null;
    String ratio = null;
    String ratioC= null;
    String depth = null;
    String support = null;

    public Entrance (String[] args) {
        this.createOptions();
        this.retrieveParameters (args);
    }

    @Override
    public void createOptions() {
        options = new Options();
        options.addOption("a", true, "App. e.g. -a Duplication");
        options.addOption("i", true, "-i /User/bin/");
        options.addOption("o", true, "-o /User/bin/");
        options.addOption("l", true, "-l 50");
        options.addOption("q", true, "-q 60");
        options.addOption("p", true, "-p /User/xujun/samtools");
        options.addOption("f", true, "-f /User/bin/");
        options.addOption("t", true, "-t 16");
        options.addOption("chr", true, "-chr 42");
        options.addOption("range", true, "-range 10");
        options.addOption("ratio", true, "-ratio 0.5");
        options.addOption("ratioC", true, "-ratio 0.2,0.4,0.6,0.8");
        options.addOption("d", true, "-d 20");
        options.addOption("s", true, "-s 2");
    }

    @Override
    public void retrieveParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            app = line.getOptionValue("a");
            if( line.hasOption( "i" ) ) {
                inputFile =line.getOptionValue("i");
            }
            if( line.hasOption( "o" ) ) {
                outputFile=line.getOptionValue("o");
            }
            if( line.hasOption( "l" ) ) {
                length=line.getOptionValue("l");
            }
            if( line.hasOption( "q" ) ) {
                quality=line.getOptionValue("q");
            }
            if( line.hasOption( "p" ) ) {
                path=line.getOptionValue("p");
            }
            if( line.hasOption( "f" ) ) {
                parameterFiles = line.getOptionValue("f");
            }
            if( line.hasOption( "t" ) ) {
                threads=line.getOptionValue("t");
            }
            if( line.hasOption( "chr" ) ) {
                chrNum=line.getOptionValue("chr");
            }
            if( line.hasOption( "range" ) ) {
                range=line.getOptionValue("range");
            }
            if( line.hasOption( "ratio" ) ) {
                ratio=line.getOptionValue("ratio");
            }
            if( line.hasOption( "ratioC" ) ) {
                ratioC=line.getOptionValue("ratioC");
            }
            if( line.hasOption( "d" ) ) {
                depth=line.getOptionValue("d");
            }
            if( line.hasOption( "s" ) ) {
                support=line.getOptionValue("s");
            }
        }
        catch(Exception e) {
            e.printStackTrace();
            System.exit(0);
        }
        if (app == null) {
            System.out.println("App does not exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
        if (app.equals(AppNames.Deletion.getName())) {
            String[] news = {this.inputFile, this.outputFile, this.length, this.quality,this.threads};
            new Deletion(news);
        }
        else if (app.equals(AppNames.Duplication.getName())) {
            String[] news ={this.inputFile,this.outputFile};
            new Duplication(news);
        }
        else if (app.equals(AppNames.Insertion.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.chrNum,this.parameterFiles,this.range,this.ratio,this.length,this.quality};
            new Insertion(news);
        }
        else if (app.equals(AppNames.InsertionNew.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.chrNum,this.parameterFiles,this.range,this.ratio,this.length,this.quality};
            new InsertionNew (news);
        }
        else if (app.equals(AppNames.InsertionOld.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.chrNum,this.parameterFiles,this.range,this.ratio,this.length,this.quality};
            new InsertionOld(news);
        }
        else if (app.equals(AppNames.MakeLib.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.parameterFiles,this.chrNum};
            new MakeLib(news);
        }
        else if (app.equals(AppNames.MakeInsLib.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.chrNum, this.quality,this.length,this.range};
            new MakeInsLib(news);
        }
        else if (app.equals(AppNames.FindFromSivm.getName())) {
            String[] news = {this.inputFile, this.parameterFiles,this.outputFile};
            new MakeInsLib(news);
        }
        else if (app.equals(AppNames.ChangeBamFlag.getName())) {
            String[] news = {this.inputFile, this.outputFile};
            new MakeInsLib(news);
        }
        else if (app.equals(AppNames.MateInsertion.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.chrNum,this.parameterFiles,this.range,this.ratio};
            new MateInsertion(news);
        }
        else if (app.equals(AppNames.UnmapInsertion.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.parameterFiles,this.ratio};
            new UnmapInsertion(news);
        }
        else if (app.equals(AppNames.Inversion.getName())) {
            String [] news ={this.inputFile,this.outputFile};
            new Inversion(news);
        }
        else if (app.equals(AppNames.GenotypeDel.getName())) {
            String news =this.parameterFiles;
            new GenotypeDel(news);
        }
        else if (app.equals(AppNames.GenotypeDelNew.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.parameterFiles,this.path,this.chrNum,this.depth,this.ratioC,this.support};
            new GenotypeDelNew(news);
        }
        else if (app.equals(AppNames.GenotypeDelNew1.getName())) {
            String[] news = {this.inputFile, this.outputFile,this.parameterFiles,this.path,this.chrNum,this.support};
            new GenotypeDelNew1(news);
        }
        else if (app.equals(AppNames.SpliteAndPE.getName())) {
            String [] news ={this.inputFile,this.outputFile,this.quality};
            new SpliteAndPE(news);
        }
        else if (app.equals(AppNames.ReliableSit.getName())) {
            String [] news ={this.inputFile,this.outputFile};
            new ReliableSit(news);
        }
        else if (app.equals(AppNames.FindURPE.getName())) {
            String [] news ={this.inputFile,this.outputFile,this.quality};
            new FindURPE(news);
        }
        else {
            System.out.println("App does not exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
    }

    @Override
    public void printIntroductionAndUsage() {
        System.out.println("Incorrect options input. Program stops.");
        System.out.println(introduction);
    }

    @Override
    public String createIntroduction() {
        StringBuilder sb = new StringBuilder();
        sb.append("\nBioinformatic toolkits of RNA-seq data is designed to simplify its usage.\n");
        sb.append("It uses two options to run its apps. \"-a\" is used to select an app. \"-f\" is used to provide a parameter file of an app.\n");
        sb.append("e.g. The command line usage of the app SiPAS-tools is: ");
        sb.append("java -Xmx100g -jar SiPAS-tools.jar -a Duplication -f parameter_parsing.txt > log.txt &\n");
        sb.append("\nAvailable apps in SiPAS-tools include,\n");
        for (int i = 0; i < AppNames.values().length; i++) {
            sb.append(AppNames.values()[i].getName()).append("\n");
        }
        sb.append("\nPlease visit https://github.com/PlantGeneticsLab/SiPAS-tools for details.\n");
        return sb.toString();
    }

    public static void main (String[] args) {new Entrance(args);}

}











