package app.Insertion;

import jdk.nashorn.internal.runtime.regexp.joni.ast.StringNode;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.File;

/**
 * @author xujun on 2022-07-01-8:58 PM
 * @project SVCaller
 */
public class Alignment {
    public void alig(){
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        SequencePair<DNASequence, NucleotideCompound> psa = null;
        DNASequence query = null; DNASequence hit = null;
        try{
            BufferedReader brC = IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/cen1.fa").getAbsolutePath());
            BufferedReader br = IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/CEN180.fa").getAbsolutePath());
            String temp = null; StringBuilder sb = new StringBuilder();
            sb.setLength(0);
            while ((temp = br.readLine()) != null){
                if (!temp.startsWith(">")){
                    sb.append(temp);
                }
            }
            query = new DNASequence(sb.toString());
            sb.setLength(0);
            while ((temp = brC.readLine()) != null){
                if (!temp.startsWith(">")){
                    sb.append(temp);
                }
            }
            hit = new DNASequence(sb.toString());
            psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);

        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


}
