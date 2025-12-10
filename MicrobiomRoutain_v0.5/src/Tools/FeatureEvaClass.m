function FeatureEvaClass()
[trainedClassifier, validationAccuracy,validationPredictions, validationScores] = trainClassifier(T_train,'Logistic',n);
[X1,Y1,T1,AUC1] = perfcurve(T_train.label,validationScores(:,1),'UIA');
display(['Train: ACC= ' num2str(validationAccuracy) ', AUC= ' num2str(AUC1)]);
figure,
plot(X1,Y1)
xlabel('1-Specificity')
ylabel('Sensitivity')
% title('ROC for Classification by Logistic (Train)')
end