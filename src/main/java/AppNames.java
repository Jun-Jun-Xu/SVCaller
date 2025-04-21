/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * Available apps in SiPAS-tools
 */
public enum AppNames implements Comparable<AppNames> {

    /**
     * Call deletion.
     */
    Deletion("Deletion"),

    /**
     * Call duplication.
     */
    Duplication("Duplication"),

    /**
     * Call Insertion.
     */
    Insertion("Insertion"),

    /**
     * Call Insertion.
     */
    InsertionNew("InsertionNew"),

    /**
     * Call Insertion.
     */
    InsertionOld("InsertionOld"),

    /**
     * Make library of Insertion.
     */
    MakeLib("MakeLib"),

    /**
     * Make library of Insertion.
     */
    MakeInsLib("MakeInsLib"),

    /**
     * Find Insertion from result of svim.
     */
    FindFromSivm("FindFromSivm"),

    /**
     * Find Insertion from result of svim.
     */
    ChangeBamFlag("ChangeBamFlag"),

    /**
     * Find insertion in unmapped reads(Both end).
     */
    MateInsertion("MateInsertion"),

    /**
     * Find insertion in unmapped reads(Both end).
     */
    UnmapInsertion("UnmapInsertion"),

    /**
     * Call inversion.
     */
    Inversion("Inversion"),

    /**
     * Genotype deletion.
     */
    GenotypeDel("GenotypeDel"),

    /**
     * Genotype deletion.
     */
    GenotypeDelNew("GenotypeDelNew"),

    /**
     * Genotype deletion.
     */
    GenotypeDelNew1("GenotypeDelNew1"),

    /**
     * Find breaking point from splite reads and unreasonble PE reads.
     */
    SpliteAndPE("SpliteAndPE"),

    /**
     * Reliable site.
     */
    ReliableSit("ReliableSit"),

    /**
     * Output unreasonble PE reads.
     */
    FindURPE("FindURPE");

    public final String name;

    AppNames(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }
}
