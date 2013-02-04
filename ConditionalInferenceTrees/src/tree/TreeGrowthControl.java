package tree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class TreeGrowthControl
{

	/**
	 * Indicates whether classification or regression is to be used.
	 */
	public boolean isClassificationUsed = true;

	/**
	 * The value of the test statistic (1 - p-value) that must be exceeded in order to implement a split.
	 * 
	 * To alter this just change the x portion of Math.log(x) to your desired value.
	 */
	public double minCriterion = Math.log(0.95);
//	public double minCriterion = -Double.MAX_VALUE;  // To grow maximal depth trees (largest negative double conceptually -Inf).

	/**
	 * The minimum sum of the weights in a node in order for it to be considered for splitting.
	 */
	public int minSplit = 20;

	/**
	 * The minimum sum of weights in a terminal node.
	 */
	public int minBucket = 7;

	/**
	 * The proportion of observations needed to establish a terminal node.
	 */
	public double minProb = 0.01;

	/**
	 * The number of input variables randomly sampled as candidates at each node.
	 */
	public int mtry = Integer.MAX_VALUE;

	/**
	 * The maximum depth that the tree can be grown to. Takes priority when determining whether to grow the tree further.
	 */
	public int maxDepth = Integer.MAX_VALUE;

	/**
	 * A boolean indicating if multiway splits for all categorical variable levels are implemented for unordered categorical variables.
	 */
	public boolean multiWay = false;

	/**
	 * The number of variables that are inspected for admissible splits if the best split doesn't meet the sample size constraints.
	 */
	public int splitTry = 2;

	/**
	 * The number of trees to grow (only use for growing a forest).
	 */
	public int numberOfTreesToGrow = 500;

	/**
	 * The variables to exclude from the growing of the tree.
	 */
	public List<String> variablesToIgnore = new ArrayList<String>();

	/**
	 * Whether or not to use replacement when selecting the examples to grow the tree (only use for growing a forest).
	 */
	public boolean isReplacementUsed = false;

	/**
	 * The fraction of observations to select if replacement is not used (only use for growing a forest).
	 */
	public double selectionFraction = 0.632;

	/**
	 * Whether to use the Bonferroni correction for the p values, or to use the univariate test type.
	 */
	public boolean isBonferronniUsed = false;

	/**
	 * Indicates whether or not a stump should be produced (a tree with three nodes).
	 */
	public boolean isStump = false;


	public TreeGrowthControl()
	{
	}

	public TreeGrowthControl(String location)
	{
		try (BufferedReader reader = Files.newBufferedReader(Paths.get(location), StandardCharsets.UTF_8))
		{
			String line = reader.readLine();
			line = line.replaceAll("\n", "");
			String ctrlVariables[] = line.split("\t");

			this.isClassificationUsed = Boolean.parseBoolean(ctrlVariables[0]);
			this.minCriterion = Double.parseDouble(ctrlVariables[1]);
			this.minSplit = Integer.parseInt(ctrlVariables[2]);
			this.minBucket = Integer.parseInt(ctrlVariables[3]);
			this.minProb = Double.parseDouble(ctrlVariables[4]);
			this.mtry = Integer.parseInt(ctrlVariables[5]);
			this.maxDepth = Integer.parseInt(ctrlVariables[6]);
			this.multiWay = Boolean.parseBoolean(ctrlVariables[7]);
			this.splitTry = Integer.parseInt(ctrlVariables[8]);
			this.numberOfTreesToGrow = Integer.parseInt(ctrlVariables[9]);
			this.variablesToIgnore = new ArrayList<String>();
			if (!ctrlVariables[10].equals(""))
			{
				String covToIgnore[] = ctrlVariables[10].split(",");
				for (String s : covToIgnore)
				{
					this.variablesToIgnore.add(s);
				}
			}
			this.isReplacementUsed = Boolean.parseBoolean(ctrlVariables[11]);
			this.selectionFraction = Double.parseDouble(ctrlVariables[12]);
			this.isBonferronniUsed = Boolean.parseBoolean(ctrlVariables[13]);
			this.isStump = Boolean.parseBoolean(ctrlVariables[14]);
		}
		catch (Exception e)
		{
			System.err.println("Error: " + e.getMessage());
			System.exit(0);
		}
	}

	/**
	 * Method to save the control object to a file.
	 * 
	 * @param location The location where the file should be saved.
	 */
	public void save(String location)
	{
		try
		{
			FileWriter outputFile = new FileWriter(location);
			BufferedWriter outputWriter = new BufferedWriter(outputFile);
			outputWriter.write(Boolean.toString(this.isClassificationUsed) + "\t");
			outputWriter.write(Double.toString(this.minCriterion) + "\t");
			outputWriter.write(Integer.toString(this.minSplit) + "\t");
			outputWriter.write(Integer.toString(this.minBucket) + "\t");
			outputWriter.write(Double.toString(this.minProb) + "\t");
			outputWriter.write(Integer.toString(this.mtry) + "\t");
			outputWriter.write(Integer.toString(this.maxDepth) + "\t");
			outputWriter.write(Boolean.toString(this.multiWay) + "\t");
			outputWriter.write(Integer.toString(this.splitTry) + "\t");
			outputWriter.write(Integer.toString(this.numberOfTreesToGrow) + "\t");
			if (!this.variablesToIgnore.isEmpty())
			{
				outputWriter.write(this.variablesToIgnore.get(0));
				for (int i = 1; i < this.variablesToIgnore.size(); i++)
				{
					outputWriter.write("," + this.variablesToIgnore.get(i));
				}
			}
			outputWriter.write("\t");
			outputWriter.write(Boolean.toString(this.isReplacementUsed) + "\t");
			outputWriter.write(Double.toString(this.selectionFraction) + "\t");
			outputWriter.write(Boolean.toString(this.isBonferronniUsed) + "\t");
			outputWriter.write(Boolean.toString(this.isStump));
			outputWriter.close();
		}
		catch (Exception e)
		{
			System.err.println(e.getStackTrace());
			System.exit(0);
		}
	}

}
