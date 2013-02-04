package tree;

import java.util.Map;


public class TestDriver
{

	public static void main(String[] args)
	{
		TreeGrowthControl ctrl = new TreeGrowthControl();
		ctrl.minCriterion = -Double.MAX_VALUE;
		ctrl.isClassificationUsed = true;
		ctrl.numberOfTreesToGrow = 100;
		Forest forest = new Forest(args[0], ctrl);
		System.out.format("The OOB error estimate is : %f\n", forest.oobErrorEstimate);

//		ProcessDataForGrowing predData = new ProcessDataForGrowing(args[0], ctrl);
//		System.out.println(forest.predict(predData, ctrl.isProbabilisticPrediction).first);
//		System.out.println(forest.predict(predData, ctrl.isProbabilisticPrediction).second);
//		for (ConditionalInferenceTree f : forest.forest)
//		{
//			Map<Integer, Object> predRes = f.predict(predData, ctrl.isProbabilisticPrediction);
//			double numberPredictions = 0;
//			double correctPRedictions = 0;
//			for (Integer i : predRes.keySet())
//			{
//				Map<String, Double> predictedClasses = (Map<String, Double>) predRes.get(i);
//				double maxProb = 0.0;
//				String maxClass = null;
//				for (String s : predictedClasses.keySet())
//				{
//					if (predictedClasses.get(s) > maxProb)
//					{
//						maxProb = predictedClasses.get(s);
//						maxClass = s;
//					}
//				}
//				numberPredictions += 1;
//				if (predData.responseData.get(maxClass).get(i) == 1)
//				{
//					correctPRedictions += 1;
//				}
//			}
//			System.out.println(correctPRedictions / numberPredictions);
//		}
//		VariableImportance varImpCalculator = new VariableImportance();
//		Map<String, Double> varImp = varImpCalculator.conditionalVariableImportance(forest, ctrl, 0.2, 2, false);
//		System.out.println(varImp.entrySet());
//		ConditionalInferenceTree ct = new ConditionalInferenceTree(args[0], ctrl);
//		ProcessDataForGrowing predData = new ProcessDataForGrowing(args[0], ctrl);
//		ct.display();
	}

}
