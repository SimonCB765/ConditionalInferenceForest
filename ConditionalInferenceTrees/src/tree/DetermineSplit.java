/**
 * 
 */
package tree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * @author Simon Bull
 *
 */
public class DetermineSplit
{

	Map<String, List<Object>> covariableData;
	Map<String, List<Double>> classData;
	int numberObservations;
	double[] weights;
	List<String> variablesToUse;
	List<String> classOrdering = new ArrayList<String>();
	int numberOfClasses;
	Map<String, String> variableTypeMapping;
	Map<String, Integer> categoricalVariableLevels;
	TreeGrowthControl ctrl;
	

	public DetermineSplit(Map<String, List<Object>> covariableData,
			Map<String, List<Double>> classData, int numberObservations, double[] weights,
			List<String> variablesToUse, Map<String, String> variableTypeMapping,
			Map<String, Integer> categoricalVariableLevels, TreeGrowthControl ctrl)
	{
		this.covariableData = covariableData;
		this.classData = classData;
		this.numberObservations = numberObservations;
		this.weights = weights;
		this.variablesToUse = variablesToUse;
		this.classOrdering.addAll(classData.keySet());
		Collections.sort(this.classOrdering);
		this.numberOfClasses = classData.keySet().size();
		this.variableTypeMapping = variableTypeMapping;
		this.categoricalVariableLevels = categoricalVariableLevels;
		this.ctrl = ctrl;
	}

	ImmutableTwoValues<Boolean, double[]> findSplitCategorical(String covariableToCheck, int minBucket,
			double minProb)
	{

		double summedWeights = 0;
		int minBucketAltered = minBucket;
		for (double i : weights)
		{
			summedWeights += i;
		}
		double summedWeightProbability = Math.ceil(summedWeights * minProb);
		if (minBucket < summedWeightProbability)
		{
			minBucketAltered = (int) summedWeightProbability;
		}
		int numberCategories = categoricalVariableLevels.get(covariableToCheck);  // The number of levels of the variable.
		// A 0 in the ith covariableSplit entry indicates that the ith category goes down the left branch.
		// A 1 in the ith covariableSplit entry indicates that the ith category goes down the right branch.
		double[] covariableSplit = new double[numberCategories];
		boolean isSplitFound = false;

		if (numberCategories == 2)
		{
			// If there are only two categories for the variable.
			int sumCategoryOne = 0;
			for (int i = 0; i < numberObservations; i++)
			{
				double observationWeight = weights[i];
				if (observationWeight == 0)
				{
					continue;
				}
				if ((int) covariableData.get(covariableToCheck).get(i) == 1)
				{
					// If the variable is from category 1 (as opposed to 2).
					sumCategoryOne += observationWeight;
				}
			}
			// Check the sample size constraints.
			// If the sum of all category 1 weights is greater than (summedWeights - minBucketAltered), this means
			// that there is not enough 'weight' in the category 2 observation to make a split.
			// For example, summedWeights = 10, minBucketaltered = 5 and sumCategoryOne = 6,
			// then the maximum 'weight' in a split made here would have 4 'weight' for the category 2 child.
			// This is less than the minimum allowed 5 weight in a terminal node as indicated by minBucketAltered.
			boolean isNotEnoughWeightCat2 = sumCategoryOne > (summedWeights - minBucketAltered);
			// If the sum of all weights for category 1 is less than minBucketAltered, then there is not enough weight
			// in the category 1 observations to make a split.
			boolean isNotEnoughWeightCat1 = sumCategoryOne < minBucketAltered;
			if (!(isNotEnoughWeightCat2 || isNotEnoughWeightCat1))
			{
				// Sample size constraints passed.
				covariableSplit[0] = 0;
				covariableSplit[1] = 1;
			}
		}
		else
		{
			// There are more than two categories.
			// Check sample size constraint (that there is enough weight using the current split).
			if (summedWeights >= minBucketAltered)
			{
				Map<String, Double> inflFuncCondExpectation = influenceFunctionCondExpectation(numberOfClasses, summedWeights);
				Map<String, Map<String, Double>> inflFuncCondCovariance = influenceFunctionCondCovariance(numberOfClasses, summedWeights, inflFuncCondExpectation);
				// A 0 in the ith entry indicates that the ith category goes down the left branch.
				// A 1 in the ith entry indicates that the ith category goes down the right branch.
				double[] possibleSplit = new double[numberCategories];
				double lastStatistic = 0.0;
				int possibleBinarySplits = 1;
				for (int i = 0; i < numberCategories - 1; i++)
				{
					// By doing i < numberCategories - 1 you get the number of binary splits as 2 ^ (numberCategories - 1).
					// This is what you want as you are fixing the first category to go down the left hand branch.
					// By doing this you cut the number of splits investigated by half, and don't run the risk
					// of missing and potential splits.
					possibleBinarySplits *= 2;
				}

				for (int i = 1; i <= possibleBinarySplits; i++)
				{
					// Start with j == 1, as you are forcing the first category to always go down the left branch,
					// and j == 0 means everything goes down the left branch.
					// This can be done as you would get the same effective splits if you allowed it to go down the right,
					// just the branch would change. For example with categories A, B and C:
					// Left : A		A,B		A,C
					// Right: B,C	C		B
					// then if you allow A to be in the right branch as well you get:
					// Left : A		A,B		A,C	||	B		B,C		C
					// Right: B,C	C		B	||	A,C		A		A,B
					// You simply have a mirror image around the ||.
					int splitter = i;
					for (int j = 1; j < numberCategories; j++)
					{
						possibleSplit[j] = splitter % 2;
						splitter /= 2;
					}

					// Compute linear statistic and the sum of weights for the 1 branch.
					// I.e. compute the stat and sum the weights using only response values which will go down
					// the 1 (right) branch if the current (possibeSplit) split is used.
					double[] linearStatValues = new double[numberOfClasses];
					double[] weightDownRightBranch = new double[1];
					weightDownRightBranch[0] = 0.0;
					for (int j = 0; j < numberObservations; j++)
					{
						// We do:
						// possibleSplit[(int) covariableData.get(covariableToCheck).get(j) - 1] == 1
						// rather than:
						// possibleSplit[(int) covariableData.get(covariableToCheck).get(j)] == 1
						// as the categories are all ordered from 1:numberOfCategories (e.g. 1,2,3 not 0,1,2).
						// This means that if the value of the covariable is 3 it corresponds to split index 2.
						if (weights[j] > 0 && possibleSplit[(int) covariableData.get(covariableToCheck).get(j) - 1] == 1)
						{
							// If the value of the covariable in this observation indicates that the observation should go down the
							// 1 (right) branch.
							weightDownRightBranch[0] += weights[j];
							for (int k = 0; k < numberOfClasses; k++)
							{
								String currentClass = classOrdering.get(k);
								linearStatValues[k] += (weights[j] * classData.get(currentClass).get(j));
							}
						}
					}

					// Check sample size constraints.
					// If the sum of all right branch weights is greater than (summedWeights - minBucketAltered), this means
					// that there is not enough 'weight' in the left branch observations to make a split.
					// For example, summedWeights = 10, minBucketAltered = 5 and weightDownRightBranch = 6,
					// then the maximum 'weight' in a split made here would have 4 'weight' for the left child.
					// This is less than the minimum allowed 5 weight in a terminal node as indicated by minBucketAltered.
					boolean isNotEnoughWeightCat2 = weightDownRightBranch[0] > (summedWeights - minBucketAltered);
					// If the sum of all weights for the right branch is less than minBucketAltered, then there is not enough weight
					// in the right branch observations to make a split.
					boolean isNotEnoughWeightCat1 = weightDownRightBranch[0] < minBucketAltered;
					if (isNotEnoughWeightCat2 || isNotEnoughWeightCat1)
					{
						continue;
					}

					double[] inflFuncCondExp = new double[numberOfClasses];
					double[] inflFuncCondCov = new double[numberOfClasses * numberOfClasses];
					int classIndex = 0;
					for (String j : classOrdering)
					{
						inflFuncCondExp[classIndex] = inflFuncCondExpectation.get(j);
						int innerClassIndex = 0;
						for (String k : classOrdering)
						{
							inflFuncCondCov[classIndex * numberOfClasses + innerClassIndex] = inflFuncCondCovariance.get(j).get(k);
							innerClassIndex++;
						}
						classIndex++;
					}
					ImmutableTwoValues<double[], double[]> condExpCov = linearStatisticCondExpCov(weightDownRightBranch,
							weightDownRightBranch, 1, numberOfClasses, summedWeights, inflFuncCondExp, inflFuncCondCov);

					// Calculate max(abs(linStatistic - explinstat) / sqrt(diag(covlinstat)))
					double[] condExpLinStat = condExpCov.first;
					double[] condCovLinStat = condExpCov.second;
					double maxAbsTestStat = maxAbsoluteTestStatistic(linearStatValues, condExpLinStat, condCovLinStat,
							numberOfClasses, 0.0);

					if (maxAbsTestStat > lastStatistic)
					{
						isSplitFound = true;
						lastStatistic = maxAbsTestStat;
						for (int j = 0; j < numberCategories; j++)
						{
							covariableSplit[j] = possibleSplit[j];
						}
					}
				}
			}
		}

		return new ImmutableTwoValues<Boolean, double[]>(isSplitFound, covariableSplit);

	}

	/**
	 * @param covariableToCheck
	 * @param minBucket
	 * @param minProb
	 * @return
	 */
	ImmutableTwoValues<Boolean, double[]> findSplitNumeric(String covariableToCheck, int minBucket, double minProb)
	{

		// Initialise the necessary variables.
		double[] splitPointFound = new double[1];
		double summedWeights = 0;
		int minBucketAltered = minBucket;
		for (double i : weights)
		{
			summedWeights += i;
		}
		double summedWeightProbability = Math.ceil(summedWeights * minProb);
		if (minBucket < summedWeightProbability)
		{
			minBucketAltered = (int) summedWeightProbability;
		}

		// Determine the order of the covariable values for each observation, along with the max value of the covariable.
		List<IndexedDoubleData> sortedCovarData = new ArrayList<IndexedDoubleData>();
		for (int i = 0; i < numberObservations; i++)
		{
			sortedCovarData.add(new IndexedDoubleData((double) covariableData.get(covariableToCheck).get(i), i));
		}
		Collections.sort(sortedCovarData);
		double maximumCovariableValue = sortedCovarData.get(numberObservations - 1).getData();

		Map<String, Double> condExpInfl = influenceFunctionCondExpectation(numberOfClasses, summedWeights);
		Map<String, Map<String, Double>> condCovInfl = influenceFunctionCondCovariance(numberOfClasses, summedWeights, condExpInfl);

		if (summedWeights >= minBucketAltered)
		{
			// The sample size constraint is passed (i.e. there are enough examples to make a split).
			// Go through the values for the covariables in order from smallest to largest.
			double[] linStatistic = new double[numberOfClasses];
			double sumOfWeightsToThisX = 0.0;  // The sum of weights for covariable values ordered <= to this one.
			double potentialSplitPoint = 0.0;  // Will be used to store a potential value of the break point while the linear statistic is being calculated.
			double bestStatValue = 0.0;
			boolean isSplitFound = false;
			for (int i = 0; i < (numberObservations - 1); i++)
			{
				double currentXValue = sortedCovarData.get(i).getData();
				int currentXRank = sortedCovarData.get(i).getIndex();
				double currentXWeight = weights[currentXRank];
				if (currentXWeight != 0)
				{
					// If the weight for the ith observation is not 0, then use it to look for a split.
					sumOfWeightsToThisX += currentXWeight;

					// Linear statistic for covariable ranks <= this rank.
					for (int j = 0; j < numberOfClasses; j++)
					{
						linStatistic[j] += weights[currentXRank] * classData.get(classOrdering.get(j)).get(currentXRank);
					}

					// Check more sample size constraints.
					if (sumOfWeightsToThisX > (summedWeights - minBucketAltered) || currentXValue >= maximumCovariableValue)
					{
						// If the sum of all weights to this point is greater than (summedWeights - minBucketAltered), this means
						// that there is not enough 'weight' remaining in the rest of the observations to make a split for any of
						// the remaining covariable values. For example, summedWeights = 10, minBucketaltered = 5 and sumOfWeightsToThisX = 6,
						// then the maximum 'weight' in a split made here would have 4 'weight' in one of the children.
						// This is less than the minimum allowed 5 weight in a terminal node as indicated by minBucketAltered.
						// If the current covariable value is >= the max value for the covariable then it is not possible to make
						// a split, as all values of the covariable would be going down one branch of any split induced.
						break;
					}
					if (sumOfWeightsToThisX < minBucketAltered)
					{
						// If this is true then the 'weight' of the obserations seen is not great enough to make a split yet.
						// You need at least minBucketAltered weight in a terminal node, so you can not split he node until the
						// observations you have gone through have reached at least that much weight.
						continue;
					}

					// Check to see that there are no values of the covariable with a rank lower than this one
					// that have the same value. For example, if currentXValue == 5 and currentXRank == 4, then
					// ensure that no covariable values with rank > 4 are equal to 5
					// (i.e. for all datapoints x, if rank(x) > 4 then value(x) > 5). If this condition is not true,
					// then loop until it is (i.e. if the values of the covariable are 1 1 1 1 2 2 3 3 etc. loop until
					//                                                                       ^   ^   ^
					int tieCheckIndex = 0;
					for (int j = i + 1; j < (numberObservations - 1); j++)
					{
						tieCheckIndex = j;
						if (weights[sortedCovarData.get(j).getIndex()] > 0)
						{
							break;
						}
					}
					if (currentXValue == sortedCovarData.get(tieCheckIndex).getData())
					{
						continue;
					}
					else
					{
						potentialSplitPoint = currentXValue;
					}

					// Calculate the conditional expectation and covariance of the linear statistic.
					double[] swxAndSwx2 = new double[1];
					swxAndSwx2[0] = sumOfWeightsToThisX;
					double[] inflFuncCondExp = new double[numberOfClasses];
					double[] inflFuncCondCov = new double[numberOfClasses * numberOfClasses];
					int classIndex = 0;
					for (String j : classOrdering)
					{
						inflFuncCondExp[classIndex] = condExpInfl.get(j);
						int innerClassIndex = 0;
						for (String k : classOrdering)
						{
							inflFuncCondCov[classIndex * numberOfClasses + innerClassIndex] = condCovInfl.get(j).get(k);
							innerClassIndex++;
						}
						classIndex++;
					}
					ImmutableTwoValues<double[], double[]> condExpCov = linearStatisticCondExpCov(swxAndSwx2,
							swxAndSwx2, 1, numberOfClasses, summedWeights, inflFuncCondExp, inflFuncCondCov);

					// Calculate max(abs(linStatistic - explinstat) / sqrt(diag(covlinstat)))
					double[] condExpLinStat = condExpCov.first;
					double[] condCovLinStat = condExpCov.second;
					double maxAbsTestStat = maxAbsoluteTestStatistic(linStatistic, condExpLinStat, condCovLinStat,
							numberOfClasses, 0.0);

					// Determine if a new best split point has been found.
					if (maxAbsTestStat > bestStatValue)
					{
						bestStatValue = maxAbsTestStat;
						splitPointFound[0] = potentialSplitPoint;
						isSplitFound = true;
					}
				}
			}
			return new ImmutableTwoValues<Boolean, double[]>(isSplitFound, splitPointFound);
		}

		return new ImmutableTwoValues<Boolean, double[]>(false, splitPointFound);

	}

	/**
	 * Calculates the conditional expectation of the influence function.
	 * 
	 * @param q
	 * @param weightSum
	 * @return
	 */
	Map<String, Double> influenceFunctionCondExpectation(int q, double weightSum)
	{

		Map<String, Double> returnValue = new HashMap<String, Double>();

		for (String key : classData.keySet())
		{
			double classWeights = 0.0;
			for (int i = 0; i < numberObservations; i++)
			{
				double classValue = classData.get(key).get(i);
				classWeights += weights[i] * classValue;
			}
			returnValue.put(key, classWeights / weightSum);
		}

		return returnValue;

	}

	/**
	 * Calculates the conditional covariance of he influence function.
	 * 
	 * @param q
	 * @param weightSum
	 * @param inflFuncCondExpectation
	 * @return
	 */
	Map<String, Map<String, Double>> influenceFunctionCondCovariance(int q, double weightSum, Map<String, Double> inflFuncCondExpectation)
	{

		Map<String, Map<String, Double>> returnValue = new HashMap<String, Map<String, Double>>();

		for (String class1 : classData.keySet())
		{
			returnValue.put(class1, new HashMap<String, Double>());
			double class1Expectation = inflFuncCondExpectation.get(class1);
			for (String class2 : classData.keySet())
			{
				double class2Expectation = inflFuncCondExpectation.get(class2);
				double currentCovariance = 0.0;
				for (int i = 0; i < numberObservations; i++)
				{
					double class1Value = classData.get(class1).get(i);
					double weightedClass1Value = weights[i] * (class1Value - class1Expectation);
					double class2Value = classData.get(class2).get(i);
					currentCovariance += weightedClass1Value * (class2Value - class2Expectation);
				}
				returnValue.get(class1).put(class2, currentCovariance / weightSum);
			}
		}

		return returnValue;

	}

	/**
	 * @param A One of the matrices.
	 * @param m # rows in A.
	 * @param n # columns in A.
	 * @param B The other matrix
	 * @param p # rows in B.
	 * @param q # columns in B.
	 * @return
	 */
	double[] kroneckerProduct(double[] A, int m, int n, double[] B, int p, int q)
	{

		int mp = m * p;  // The # rows in the result matrix.
		int ip;  // Used to index the first element of each row in the result matrix.
		int jq;  // Used to index the first element of each column in the result matrix.
		double[] returnValue = new double[mp * n * q];  // Initialise the mp by nq result matrix.
		double y;

	    for (int i = 0; i < m; i++)
	    {
	    	ip = i * p;
	    	for (int j = 0; j < n; j++)
	    	{
	    		jq = j * q;
	    		y = A[j * m + i];
	    		for (int k = 0; k < p; k++)
	    		{
	    			for (int l = 0; l < q; l++)
	    			{
	    				returnValue[(jq + l) * mp + ip + k] = y * B[l * p + k];
	                }
	            }
	        }
	    }

		return returnValue;

	}

	/**
	 * @param covariableData
	 * @param p
	 * @return
	 */
	Map<String, Map<Integer, Double>> linearStatisticCategorical(List<Object> covariableData, int p)
	{

		Map<String, Map<Integer, Double>> returnValue = new HashMap<String, Map<Integer, Double>>();

		for (String key : classData.keySet())
		{
			Map<Integer, Double> classWeights = new HashMap<Integer, Double>();
			for (int i = 1; i <= p; i++)
			{
				classWeights.put(i, 0.0);
			}
			for (int i = 0; i < numberObservations; i++)
			{
				Integer observedValue = (Integer) covariableData.get(i);
				double classValue = classData.get(key).get(i);
				double weightedValue = classValue * weights[i];
				classWeights.put(observedValue, classWeights.get(observedValue) + weightedValue);
			}
			returnValue.put(key, classWeights);
		}

		return returnValue;

	}

	/**
	 * Calculates the linear statistic.
	 * 
	 * For a numeric covariable this is essentially the weighted sum of the observations that are
	 * classified as a given class. For example, you have:
	 * covariable	class1 class2
	 * 2			1	0
	 * 3			0	1
	 * 5			0	1
	 * 7			1	0
	 * 4			0	1
	 * Then the value for class1 is 9 and the value for class2 is 12.
	 * The Map returned will be {'class1' : 9, 'class2' : 12}.
	 * 
	 * @param covariableData
	 * @param q
	 * @return
	 */
	Map<String, Double> linearStatisticNumeric(List<Object> covariableData)
	{

		Map<String, Double> returnValue = new HashMap<String, Double>();

		for (String key : classData.keySet())
		{
			double classWeights = 0.0;
			for (int i = 0; i < numberObservations; i++)
			{
				Double observedValue = (Double) covariableData.get(i);
				Double weightedValue = observedValue * weights[i];
				double classValue = classData.get(key).get(i);
				classWeights += weightedValue * classValue;
			}
			returnValue.put(key, classWeights);
		}

		return returnValue;

	}

	/**
	 * @param swx
	 * @param swx2
	 * @param p
	 * @param q
	 * @param weightSum
	 * @param inflFuncCondExpectation
	 * @param influenceFuncCondCovariance
	 * @return
	 */
	ImmutableTwoValues<double[], double[]> linearStatisticCondExpCov(double[] swx, double[] swx2, int p, int q, double weightSum,
			double[] inflFuncCondExpectation, double[] influenceFuncCondCovariance)
	{

		int pq = p * q;

		// Calculate the expectation of the linear statistic.
		double[] expLinearStat = new double[pq];
		for (int i = 0; i < p; i++)
		{
			for (int j = 0; j < q; j++)
			{
				expLinearStat[j * p + i] = swx[i] * inflFuncCondExpectation[j];
			}
		}

		// Calculate the conditional covariance of the linear statistic.
		double sumOfWeights = weightSum;
		double modifiedSumOfWeights = sumOfWeights / (sumOfWeights - 1);
		double inverseSumOfWeights = 1 / (sumOfWeights - 1);
		double[] covLinearStat = kroneckerProduct(influenceFuncCondCovariance, q, q, swx2, p, p);
		double[] covInfXSwx = kroneckerProduct(influenceFuncCondCovariance, q, q, swx, p, 1);
		double[] tempCov = kroneckerProduct(covInfXSwx, pq, q, swx, 1, p);
		for (int i = 0; i < (pq * pq); i++)
		{
			covLinearStat[i] = (modifiedSumOfWeights * covLinearStat[i]) - (inverseSumOfWeights * tempCov[i]);
		}

		return new ImmutableTwoValues<double[], double[]>(expLinearStat, covLinearStat);

	}

	/**
	 * @param linStatistic
	 * @param condExpLinStat
	 * @param condCovLinStat
	 * @param pq
	 * @param tolerance
	 * @return
	 */
	double maxAbsoluteTestStatistic(double[] linStatistic, double[] condExpLinStat, double[] condCovLinStat,
			int pq, double tolerance)
	{
		double resultValue = 0.0;
		double tmp = 0.0;
		double sd;

		for (int i = 0; i < pq; i++)
		{
			sd = condCovLinStat[i * pq + i];
			if (sd > tolerance)
			{
				tmp = Math.abs((linStatistic[i] - condExpLinStat[i]) / Math.sqrt(sd));
			}
			if (tmp > resultValue)
			{
				resultValue = tmp;
			}
		}

		return resultValue;
	}

	/**
	 * @param linStatVals
	 * @param expLinearStat
	 * @param covLinearStat
	 * @return
	 */
	double[] quadraticPValue(double[] linStatVals, double[] expLinearStat, double[] covLinearStat)
	{

		// Create a column vector of the difference between the linear statistic values and the expected values.
		double[] linExpDifTmp = new double[linStatVals.length];
		for (int i = 0; i < linStatVals.length; i++)
		{
			linExpDifTmp[i] = linStatVals[i] - expLinearStat[i];
		}
		RealMatrix linExpDif = new Array2DRowRealMatrix(linExpDifTmp);

		// Create a 2D covariance matrix.
		int columnDimension = linStatVals.length;
		int rowDimension = covLinearStat.length / columnDimension;
		RealMatrix covarianceMatrix = new Array2DRowRealMatrix(rowDimension, columnDimension);
		for (int i = 0; i < rowDimension; i++)
		{
			for (int j = 0; j < columnDimension; j++)
			{
				covarianceMatrix.addToEntry(i, j, covLinearStat[i * columnDimension + j]);
			}
		}

		// Calculate the SVD of the covariance matrix.
		SingularValueDecomposition svd = new SingularValueDecomposition(covarianceMatrix);
		double tolerance = Math.sqrt(2.220446e-16);  // Was like this in R code : Math.sqrt(Double.MIN_VALUE);
		                                             // but instead I will use the value from R rather than the method.
		RealMatrix sigma = svd.getS();  // The sigma portion of the SVD.
		RealMatrix decompositionU = svd.getU();
		RealMatrix decompositionV = svd.getV();

		// Determine where the Sigma values from the SVD can be deemed to be positive.
		List<Boolean> positiveSigma = new ArrayList<Boolean>();
		List<Integer> tmpPosIndices = new ArrayList<Integer>();  // Records the indices where positiveSigma is true.
		for (int i = 0; i < sigma.getColumnDimension(); i++)
		{
			boolean posSig = sigma.getEntry(i, i) > Math.max(tolerance * sigma.getEntry(0, 0), 0);
			positiveSigma.add(posSig);
			if (posSig)
			{
				tmpPosIndices.add(i);
			}
		}
		int numberOfPositives = Collections.frequency(positiveSigma, true);
		int[] positiveSigmaIndices = new int[tmpPosIndices.size()];  // Used to slice the matrices.
		for (int i = 0; i < tmpPosIndices.size(); i++)
		{
			positiveSigmaIndices[i] = (int) tmpPosIndices.get(i);
		}

		RealMatrix XPlus = new Array2DRowRealMatrix();
		if (numberOfPositives == positiveSigma.size())
		{
			RealMatrix intermediary = new Array2DRowRealMatrix(rowDimension, columnDimension);
			RealMatrix transposeU = decompositionU.transpose();
			for (int i = 0; i < rowDimension; i++)
			{
				double inverseSigmaValue = 1 / sigma.getEntry(i, i);
				for (int j = 0; j < columnDimension; j++)
				{
					intermediary.setEntry(i, j, inverseSigmaValue * transposeU.getEntry(i, j));
				}
			}
			XPlus = decompositionV.multiply(intermediary);
		}
		else if (numberOfPositives == 0)
		{
			// A matrix of all 0s with the same dimensions as the covariance matrix.
			XPlus = new Array2DRowRealMatrix(rowDimension, columnDimension);
		}
		else
		{
			int[] rowsToUse = new int[decompositionU.getRowDimension()];
			for (int i = 0; i < decompositionU.getRowDimension(); i++)
			{
				rowsToUse[i] = i;
			}
			RealMatrix vSubMatrix = decompositionV.getSubMatrix(rowsToUse, positiveSigmaIndices);
			RealMatrix transposeUSubMatrix = decompositionU.getSubMatrix(rowsToUse, positiveSigmaIndices).transpose();
			RealMatrix intermediary = new Array2DRowRealMatrix(transposeUSubMatrix.getRowDimension(), transposeUSubMatrix.getColumnDimension());
			for (int i: positiveSigmaIndices)
			{
				double inverseSigmaValue = 1 / sigma.getEntry(i, i);
				for (int j = 0; j < rowDimension; j++)
				{
					// Transposed the U matrix, so raher than rowDim by colim it is now colDim by rowDim.
					intermediary.setEntry(i, j, inverseSigmaValue * transposeUSubMatrix.getEntry(i, j));
				}
			}
			XPlus = vSubMatrix.multiply(intermediary);
		}
		int degreesOfFreedom = numberOfPositives;  // The degrees of freedom for use with the chi squared test.

		// Finish calculating the p values.
		RealMatrix crossprod = (linExpDif.transpose()).multiply(XPlus);
		double testStatistic = crossprod.multiply(linExpDif).getEntry(0, 0);
		testStatistic = Math.max(testStatistic, 0.0);
		double pValue = 0.0;
		if (degreesOfFreedom == 0 && testStatistic == 0.0)
		{
			// If there are no degrees of freedom this means that the node only contains examples of one class.
			pValue = -Double.MAX_VALUE;
		}
		else
		{
			ChiSquaredDistribution chi = new ChiSquaredDistribution(degreesOfFreedom);
			pValue = Math.log(chi.cumulativeProbability(testStatistic));
			testStatistic = Math.log(testStatistic);
		}

		double[] returnValue = new double[3];
		returnValue[0] = testStatistic;
		returnValue[1] = pValue;
		returnValue[2] = pValue;
		return returnValue;

	}

	/**
	 * Ranks the covariables based on their p values.
	 * 
	 * Returns only those covariables which have a p value that is significant at the level requested. 
	 * 
	 * @return
	 */
	Map<String, double[]> rankVariablesForSplitting()
	{

		Map<String, double[]> returnMap = new HashMap<String, double[]>();
		boolean isBonferroniUsed = this.ctrl.isBonferronniUsed;

		double weightSum = 0;
		for (double i : weights)
		{
			weightSum += i;
		}

		// For every variable under consideration, calculate the variables p value.
		for (String i : variablesToUse)
		{
			int p;
			double[] swx;
			double[] swx2;
			// For categorical variables each row represents a class and each column a different category.
			// The entries hold the number of observations where category c corresponds to class k, e.g.
			// 		cat1	cat2	cat3
			// cl1	4		5		9
			// cl2	7		5		8
			// cl3	1		0		6
			double[] linStatVals;
			

			if (variableTypeMapping.get(i).equals("c"))
			{
				// If the variable is categorical.
				p = categoricalVariableLevels.get(i);  // The number of levels of the variable.
				linStatVals = new double[p * numberOfClasses];
				List<Object> observedVariableValues = covariableData.get(i);
				ImmutableTwoValues<double[], double[]> weightSums = swxCategorical(observedVariableValues, p);
				swx = weightSums.first;
				swx2 = weightSums.second;
				Map<String, Map<Integer, Double>> linearStatisticValues = linearStatisticCategorical(observedVariableValues, p);
				// Convert the linear statistic Map to a simpler array representation.
				int classIndex = 0;
				for (String j : classOrdering)
				{
					for (int k = 0; k < p; k++)
					{
						linStatVals[classIndex * p + k] = linearStatisticValues.get(j).get(k + 1);
					}
					classIndex++;
				}
			}
			else
			{
				// If the variable is numeric.
				p = 1;  // The number of levels of the variable.
				linStatVals = new double[p * numberOfClasses];
				List<Object> observedVariableValues = covariableData.get(i);
				ImmutableTwoValues<double[], double[]> weightSums = swxNumeric(observedVariableValues);
				swx = weightSums.first;
				swx2 = weightSums.second;
				Map<String, Double> linearStatisticValues = linearStatisticNumeric(observedVariableValues);
				// Convert the linear statistic Map to a simpler array representation.
				int classIndex = 0;
				for (String j : classOrdering)
				{
					linStatVals[classIndex] = linearStatisticValues.get(j);
					classIndex++;
				}
			}

			Map<String, Double> inflFuncCondExpectation = influenceFunctionCondExpectation(numberOfClasses, weightSum);
			Map<String, Map<String, Double>> influenceFuncCondCovariance = influenceFunctionCondCovariance(numberOfClasses, weightSum, inflFuncCondExpectation);

			// Convert the influence function Maps to simpler array representations of matrices for further calculations.
			double[] inflFuncCondExp = new double[p * numberOfClasses];
			double[] inflFuncCondCov = new double[(p * numberOfClasses) * (p * numberOfClasses)];
			int classIndex = 0;
			for (String j : classOrdering)
			{
				inflFuncCondExp[classIndex] = inflFuncCondExpectation.get(j);
				int innerClassIndex = 0;
				for (String k : classOrdering)
				{
					inflFuncCondCov[classIndex * numberOfClasses + innerClassIndex] = influenceFuncCondCovariance.get(j).get(k);
					innerClassIndex++;
				}
				classIndex++;
			}

			ImmutableTwoValues<double[], double[]> linStatResults = linearStatisticCondExpCov(swx, swx2, p, numberOfClasses, weightSum, inflFuncCondExp,
					inflFuncCondCov);
			double[] expLinearStat = linStatResults.first;
			double[] covLinearStat = linStatResults.second;

			// Now use the calculated values to determine the p values.
			double[] testStatAndPValue = quadraticPValue(linStatVals, expLinearStat, covLinearStat);

			if (isBonferroniUsed)
			{
				// Correct the p value using the Bonferroni correction.
				// Multiply the p value by the number of covariables (equivalent to the number of tests performed).
				testStatAndPValue[1] = testStatAndPValue[1] * variablesToUse.size();
				testStatAndPValue[2] = testStatAndPValue[2] * variablesToUse.size();
			}

			returnMap.put(i, testStatAndPValue);

		}

		return returnMap;

	}

	/**
	 * Helper function for the calculation of the conditional expectation and covariance for categorical covariables.
	 * 
	 * @param covariableData
	 * @param numberCategories
	 */
	ImmutableTwoValues<double[], double[]> swxCategorical(List<Object> covariableData, int numberCategories)
	{
	
		double[] swx = new double[numberCategories];
		double[] swx2 = new double[numberCategories * numberCategories];

		// Expectation portion.
		for (int i = 0; i < numberObservations; i++)
		{
			Integer observedValue = (Integer) covariableData.get(i);
			swx[observedValue - 1] += weights[i];
		}

		// Covariance portion.
		for (int i = 0; i < numberCategories; i++)
		{
			swx2[i * numberCategories + i] = swx[i];
		}

		return new ImmutableTwoValues<double[], double[]>(swx, swx2);
	
	}

	/**
	 * Helper function for the calculation of the conditional expectation and covariance for numeric covariables.
	 * 
	 * Calculates
	 * swx = sigma(i = 1:n) (weight[i] * observedValue[i])
	 * and
	 * swx2 = sigma(i = 1:n) (the Kronecker product of (weight[i] * observedValue[i]) and transpose(weight[i] * observedValue[i])))
	 * 
	 * @param covariableData
	 * @return
	 */
	ImmutableTwoValues<double[], double[]> swxNumeric(List<Object> covariableData)
	{

		double[] swx = new double[1];
		double[] swx2 = new double[1];

		for (int i = 0; i < numberObservations; i++)
		{
			Double observedValue = (Double) covariableData.get(i);
			swx[0] += weights[i] * observedValue;
			swx2[0] += weights[i] * observedValue * observedValue;
			
		}

		return new ImmutableTwoValues<double[], double[]>(swx, swx2);

	}

}
