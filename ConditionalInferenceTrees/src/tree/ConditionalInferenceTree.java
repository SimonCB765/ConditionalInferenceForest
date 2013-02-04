/**
 * 
 */
package tree;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author		Simon Bull
 * @version		1.0
 * @since		2012-11-30
 */
public class ConditionalInferenceTree
{

	/**
	 * The tree that has been grown.
	 */
	Node condInfTree = null;

	/**
	 * The object recording the control parameters for the tree.
	 */
	TreeGrowthControl ctrl;

	/**
	 * 
	 */
	ProcessDataForGrowing processedData;

	public ConditionalInferenceTree(String dataFileToGrowFrom)
	{
		TreeGrowthControl ctrl = new TreeGrowthControl();
		this.ctrl = ctrl;
		ProcessDataForGrowing processedData = new ProcessDataForGrowing(dataFileToGrowFrom, this.ctrl);
		this.processedData = processedData;
		controlTreeGrowth();
	}

	public ConditionalInferenceTree(String dataFileToGrowFrom, TreeGrowthControl ctrl)
	{
		this.ctrl = ctrl;
		ProcessDataForGrowing processedData = new ProcessDataForGrowing(dataFileToGrowFrom, this.ctrl);
		this.processedData = processedData;
		controlTreeGrowth();
	}

	public ConditionalInferenceTree(String dataFileToGrowFrom, Double[] weights)
	{
		TreeGrowthControl ctrl = new TreeGrowthControl();
		this.ctrl = ctrl;
		ProcessDataForGrowing processedData = new ProcessDataForGrowing(dataFileToGrowFrom, this.ctrl);
		this.processedData = processedData;
		controlTreeGrowth(weights);
	}

	public ConditionalInferenceTree(String dataFileToGrowFrom, TreeGrowthControl ctrl, Double[] weights)
	{
		this.ctrl = ctrl;
		ProcessDataForGrowing processedData = new ProcessDataForGrowing(dataFileToGrowFrom, this.ctrl);
		this.processedData = processedData;
		controlTreeGrowth(weights);
	}

	public ConditionalInferenceTree(ProcessDataForGrowing processedData)
	{
		TreeGrowthControl ctrl = new TreeGrowthControl();
		this.ctrl = ctrl;
		this.processedData = processedData;
		controlTreeGrowth();
	}

	public ConditionalInferenceTree(ProcessDataForGrowing processedData, TreeGrowthControl ctrl)
	{
		this.ctrl = ctrl;
		this.processedData = processedData;
		controlTreeGrowth();
	}

	public ConditionalInferenceTree(ProcessDataForGrowing processedData, Double[] weights)
	{
		TreeGrowthControl ctrl = new TreeGrowthControl();
		this.ctrl = ctrl;
		this.processedData = processedData;
		controlTreeGrowth(weights);
	}

	public ConditionalInferenceTree(ProcessDataForGrowing processedData, TreeGrowthControl ctrl, Double[] weights)
	{
		this.ctrl = ctrl;
		this.processedData = processedData;
		controlTreeGrowth(weights);
	}

	void controlTreeGrowth()
	{
		controlTreeGrowth(null);
	}

	/**
	 * Controls the growth of the tree.
	 */
	void controlTreeGrowth(Double[] potentialWeights)
	{

		double[] weights = new double[this.processedData.numberObservations];
		if (potentialWeights == null)
		{
			// If the weight array is not supplied, then initialise the weights array to an array of all 1s.
			Arrays.fill(weights, 1.0);
		}
		else
		{
			// If weights were supplied make sure they are all positive doules, and that the array has
			// the correct length.
			if (potentialWeights.length != this.processedData.numberObservations)
			{
				System.out.format("The length of the weight array must be the same as the number of observations.\nThere are %d observations, but the supplied weight array contains %d elements.", this.processedData.numberObservations, potentialWeights.length);
				System.exit(0);
			}
			for (int i = 0; i < this.processedData.numberObservations; i++)
			{
				if (potentialWeights[i] instanceof Double)
				{
					weights[i] = (double) potentialWeights[i];
				}
				else
				{
					System.out.format("All weights must be doubles. Element number %d is not.", i);
					System.exit(0);
				}
			}
		}

		// Grow the tree.
		this.condInfTree = growTree(this.processedData.covariableData, this.processedData.responseData,
				this.processedData.numberObservations, weights, ctrl, 0);

	}

	/**
	 * Returns a set of the values for the split points in the tree where the node is split on the variable covariable.
	 * @param covariable
	 * @return
	 */
	List<Double> cutpoints(String covariable)
	{
		return this.condInfTree.cutpoints(covariable);
	}

	/**
	 * Displays the tree.
	 */
	public void display()
	{

		this.condInfTree.display();

	}

	/**
	 * @param covariableData
	 * @param responseData
	 * @param numberObservations
	 * @param weights
	 * @param ctrl
	 * @param currentDepth
	 * @return
	 */
	Node growTree(Map<String, List<Object>> covariableData, Map<String, List<Double>> responseData,
			int numberObservations, double[] weights, TreeGrowthControl ctrl, int currentDepth)
	{

		// Determine the counts of each class in the current node.
		int numberOfObservationsInNode = 0;
		double meanResponseValue = 0.0;
		Map<String, Integer> classCountsForNode = new HashMap<String, Integer>();
		for (Map.Entry<String, List<Double>> entry : responseData.entrySet())
		{
			String key = entry.getKey();
			List<Double> value = entry.getValue();
			int entrySum = 0;
			for (int i = 0; i < numberObservations; i++)
			{
				if (weights[i] != 0)
				{
					if (ctrl.isClassificationUsed)
					{
						entrySum += value.get(i);
					}
					else
					{
						// Determine the mean value of the response variable for the observations in this node.
						entrySum += 1;
						meanResponseValue += value.get(i);
					}
				}
			}
			classCountsForNode.put(key, entrySum);
			numberOfObservationsInNode += entrySum;
		}
		meanResponseValue /= numberOfObservationsInNode;

		//**********************************************
		// Check whether growth should stop.
		//**********************************************
		if (currentDepth == ctrl.maxDepth)
		{
			// The maximum depth permissible has been reached. A terminal node must therefore be created.
			if (ctrl.isClassificationUsed)
			{
				Map<String, Double> weightedClassCounts = new HashMap<String, Double>();
				for (String s : this.processedData.responseData.keySet())
				{
					weightedClassCounts.put(s, 0.0);
				}
				for (int i = 0; i < this.processedData.numberObservations; i++)
				{
					double obsWeight = weights[i];
					if (obsWeight == 0)
					{
						continue;
					}
					String obsClass = "";
					for (String s : this.processedData.responseData.keySet())
					{
						if (this.processedData.responseData.get(s).get(i) == 1)
						{
							obsClass = s;
						}
					}
					weightedClassCounts.put(obsClass, weightedClassCounts.get(obsClass) + obsWeight);
				}
				return new TerminalClassificationNode(classCountsForNode, currentDepth, weightedClassCounts);
			}
			else
			{
				return new TerminalRegressionNode(currentDepth, numberOfObservationsInNode, meanResponseValue);
			}
		}
		if (currentDepth > 0 && ctrl.isStump)
		{
			// If you have already created the root and a stump is to be produced, then this node must be terminal.
			if (ctrl.isClassificationUsed)
			{
				Map<String, Double> weightedClassCounts = new HashMap<String, Double>();
				for (String s : this.processedData.responseData.keySet())
				{
					weightedClassCounts.put(s, 0.0);
				}
				for (int i = 0; i < this.processedData.numberObservations; i++)
				{
					double obsWeight = weights[i];
					if (obsWeight == 0)
					{
						continue;
					}
					String obsClass = "";
					for (String s : this.processedData.responseData.keySet())
					{
						if (this.processedData.responseData.get(s).get(i) == 1)
						{
							obsClass = s;
						}
					}
					weightedClassCounts.put(obsClass, weightedClassCounts.get(obsClass) + obsWeight);
				}
				return new TerminalClassificationNode(classCountsForNode, currentDepth, weightedClassCounts);
			}
			else
			{
				return new TerminalRegressionNode(currentDepth, numberOfObservationsInNode, meanResponseValue);
			}
		}
		double weightSum = 0;
		for (double i : weights)
		{
			weightSum += i;
		}
		if (weightSum < ctrl.minSplit)
		{
			// The sum of the weights remaining is less than the sum required to perform a split.
			// A terminal node must therefore be created.
			if (ctrl.isClassificationUsed)
			{
				Map<String, Double> weightedClassCounts = new HashMap<String, Double>();
				for (String s : this.processedData.responseData.keySet())
				{
					weightedClassCounts.put(s, 0.0);
				}
				for (int i = 0; i < this.processedData.numberObservations; i++)
				{
					double obsWeight = weights[i];
					if (obsWeight == 0)
					{
						continue;
					}
					String obsClass = "";
					for (String s : this.processedData.responseData.keySet())
					{
						if (this.processedData.responseData.get(s).get(i) == 1)
						{
							obsClass = s;
						}
					}
					weightedClassCounts.put(obsClass, weightedClassCounts.get(obsClass) + obsWeight);
				}
				return new TerminalClassificationNode(classCountsForNode, currentDepth, weightedClassCounts);
			}
			else
			{
				return new TerminalRegressionNode(currentDepth, numberOfObservationsInNode, meanResponseValue);
			}
		}
		

		//**********************************************
		// Determine p values for the covariables.
		//**********************************************
		// First select the variables to check for this split.
		// Do this by generating a random shuffle of the covariables,
		// and then selecting the first mtry covariables as the ones to use.
		Set<String> covariablesAvailable = covariableData.keySet();
		List<String> shuffledCovariables = new ArrayList<String>(covariablesAvailable);
		Collections.shuffle(shuffledCovariables);
		int numVarsToSelect = Math.min(covariablesAvailable.size(), ctrl.mtry);
		List<String> variablesToUse = shuffledCovariables.subList(0, numVarsToSelect);

		// Next determine the p values for the covariables chosen to be evaluated for the split.
		DetermineSplit splitCalculator = new DetermineSplit(covariableData, responseData, numberObservations,
				weights, variablesToUse, this.processedData.variableTypeMapping, this.processedData.categoricalVariableLevels,
				this.ctrl);
		Map<String, double[]> covarPValues = splitCalculator.rankVariablesForSplitting();

		//**********************************************
		// Try to find a split.
		//**********************************************
		int variablesTried = 0;
		boolean isSplitFound = false;
		boolean isPotentiallySignificant = true;
		double[] splitValue = null;
		ImmutableTwoValues<Boolean, double[]> splitFound = null;
		String covarWithMaxPValue = "";
		// Loop while a split has not been found, there may be a significant covariable and a split should still be searched for.
		while (variablesTried < ctrl.splitTry && !isSplitFound && isPotentiallySignificant)
		{
			// Search for a covariable with the p value that is the greater than the minimum acceptable p value,
			// and is the greatest of all covariable p values.
			double maxPValuefound = -Double.MAX_VALUE;
			boolean isSignificantCovar = false;  // Indicates whether or not there is a covariable that should be split on.
			for (String i : covarPValues.keySet())
			{
				if (covarPValues.get(i)[1] > ctrl.minCriterion)
				{
					isSignificantCovar = true;
					if (covarPValues.get(i)[1] > maxPValuefound)
					{
						maxPValuefound = covarPValues.get(i)[1];
						covarWithMaxPValue = i;
					}
				}
			}

			// If there is a significant covariable, then look for a valid split in that covariable.
			// If there isn't then record this so the loop terminates.
			if (isSignificantCovar)
			{
				// Check the covariable for a valid split.
				if (this.processedData.variableTypeMapping.get(covarWithMaxPValue).equals("c"))
				{
					// If the covariable with the maximum p value is a categorical variable.
					splitFound = splitCalculator.findSplitCategorical(covarWithMaxPValue, ctrl.minBucket, ctrl.minProb);
				}
				else
				{
					// If the covariable with the maximum p value is not a categorical variable.
					splitFound = splitCalculator.findSplitNumeric(covarWithMaxPValue, ctrl.minBucket, ctrl.minProb);
				}
				if (splitFound.first)
				{
					isSplitFound = true;
					splitValue = splitFound.second;
				}
				else
				{
					variablesTried++;
					covarPValues.get(covarWithMaxPValue)[1] = -Double.MAX_VALUE;
				}
			}
			else
			{
				isPotentiallySignificant = false;
			}
		}

		// Return the Node and continue building the tree.
		if (isSplitFound)
		{
			// If a valid split was found, then generate a non-terminal node and recurse through its children.
			if (!this.processedData.variableTypeMapping.get(covarWithMaxPValue).equals("c"))
			{
				// If the variable with the max p value is not categorical.
				double[] leftWeights = new double[numberObservations];
				double[] rightWeights = new double[numberObservations];
				for (int i = 0; i < numberObservations; i++)
				{
					if ((double) covariableData.get(covarWithMaxPValue).get(i) > splitValue[0])
					{
						rightWeights[i] = weights[i];
						leftWeights[i] = 0;
					}
					else
					{
						rightWeights[i] = 0;
						leftWeights[i] = weights[i];
					}
				}
				Node leftChild = growTree(covariableData, responseData, numberObservations, leftWeights, ctrl, currentDepth + 1);
				Node rightChild = growTree(covariableData, responseData, numberObservations, rightWeights, ctrl, currentDepth + 1);
				if (ctrl.isClassificationUsed)
				{
					return new NonTerminalNumericClassificationNode(currentDepth, covarWithMaxPValue, splitValue[0], leftChild,
							rightChild, classCountsForNode, covarPValues);
				}
				else
				{
					return new NonTerminalNumericRegressionNode(currentDepth, covarWithMaxPValue, splitValue[0], leftChild,
							rightChild, covarPValues, numberOfObservationsInNode);
				}
			}
			else
			{
				// If the variable with the max p value is categorical.
				// Initialise the weight arrays for each child node.
				List<double[]> weightArrays = new ArrayList<double[]>();
				if (ctrl.multiWay)
				{
					// Multiway splits.
//					int covariableCategories = categoricalVariableLevels.get(covarWithMaxPValue);
//					for (int i = 0; i < covariableCategories; i++)
//					{
//						weightArrays.add(new double[numberObservations]);
//					}
				}
				else
				{
					// Binary split.
					weightArrays.add(new double[numberObservations]);
					weightArrays.add(new double[numberObservations]);
				}

				// Determine the weights for each child node.
				for (int i = 0; i < numberObservations; i++)
				{
					int covariableValue = (int) covariableData.get(covarWithMaxPValue).get(i);  // The value of the covariable for the observation.
					int childIndex = (int) splitValue[covariableValue - 1];  // The branch (child node) that the covariable value indicates that this observation should go down.
					for (int j = 0; j < weightArrays.size(); j++)
					{
						if (j == childIndex)
						{
							// If the branch/child index is the same as the weight array index, then this is
							// the weight array/branch that the observation should go with.
							weightArrays.get(j)[i] = weights[i];
						}
						else
						{
							weightArrays.get(j)[i] = 0;
						}
					}
				}

				// Create all the child nodes for the current node.
				List<Node> childrenOfCurrentNode = new ArrayList<Node>();
				for (double[] wa : weightArrays)
				{
					childrenOfCurrentNode.add(growTree(covariableData, responseData, numberObservations,
							wa, ctrl, currentDepth + 1));
				}

				if (ctrl.isClassificationUsed)
				{
					return new NonTerminalCategoricalClassificationNode(currentDepth, covarWithMaxPValue, splitValue,
							childrenOfCurrentNode, classCountsForNode, covarPValues);
				}
				else
				{
					return new NonTerminalCategoricalRegressionNode(currentDepth, covarWithMaxPValue, splitValue,
							childrenOfCurrentNode, covarPValues, numberOfObservationsInNode);
				}
			}
		}
		else
		{
			// If there was no valid split found, then return a terminal node.
			if (ctrl.isClassificationUsed)
			{
				Map<String, Double> weightedClassCounts = new HashMap<String, Double>();
				for (String s : this.processedData.responseData.keySet())
				{
					weightedClassCounts.put(s, 0.0);
				}
				for (int i = 0; i < this.processedData.numberObservations; i++)
				{
					double obsWeight = weights[i];
					if (obsWeight == 0)
					{
						continue;
					}
					String obsClass = "";
					for (String s : this.processedData.responseData.keySet())
					{
						if (this.processedData.responseData.get(s).get(i) == 1)
						{
							obsClass = s;
						}
					}
					weightedClassCounts.put(obsClass, weightedClassCounts.get(obsClass) + obsWeight);
				}
				return new TerminalClassificationNode(classCountsForNode, currentDepth, weightedClassCounts);
			}
			else
			{
				return new TerminalRegressionNode(currentDepth, numberOfObservationsInNode, meanResponseValue);
			}
		}

	}

	Set<String> nonTerminalCovariablesInTree()
	{
		return this.condInfTree.getChildren();
	}

	Map<Integer, Object> predict(ProcessDataForGrowing predData)
	{
		List<Integer> observationsToPredict = new ArrayList<Integer>();
		for (int i = 0; i < predData.numberObservations; i++)
		{
			observationsToPredict.add(i);
		}
		return predict(predData, observationsToPredict);
	}

	Map<Integer, Object> predict(ProcessDataForGrowing predData, List<Integer> observationsToPredict)
	{
		Map<Integer, Object> returnValue = new HashMap<Integer, Object>();
		if (this.condInfTree == null)
		{
			returnValue = null;
		}
		else
		{
			for (Integer i : observationsToPredict)
			{
				// The current observation is a mapping from the covariable names to their values.
				Map<String, Object> currentObservation = new HashMap<String, Object>();
				for (String s : predData.covariableData.keySet())
				{
					currentObservation.put(s, predData.covariableData.get(s).get(i));
				}
				returnValue.put(i, condInfTree.predict(currentObservation));
			}
		}

		return returnValue;
	}

	void save(String location)
	{
		try
		{
			FileWriter outputFile = new FileWriter(location);
			BufferedWriter outputWriter = new BufferedWriter(outputFile);
			outputWriter.close();
		}
		catch (Exception e)
		{
			System.err.println(e.getStackTrace());
			System.exit(0);
		}
	}

}
