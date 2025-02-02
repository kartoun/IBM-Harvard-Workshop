drop view v_EMRBots_Lab_Helper_1
go
Create view v_EMRBots_Lab_Helper_1
as
select t3.PatientID, t3.AdmissionID, t3.LabName, t3.Most_Recent_Value_During_Admission
from
(
select t2.*
,(select LabValue where LabDateTime = Most_Recent_Date) as Most_Recent_Value_During_Admission
from
(
select t1.*
,(case when LabName = 'CBC: ABSOLUTE NEUTROPHILS' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'CBC: ABSOLUTE NEUTROPHILS' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID) 
		when LabName = 'CBC: HEMATOCRIT' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'CBC: HEMATOCRIT' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'CBC: HEMOGLOBIN' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'CBC: HEMOGLOBIN' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'CBC: PLATELET COUNT' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'CBC: PLATELET COUNT' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'CBC: RED BLOOD CELL COUNT' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'CBC: RED BLOOD CELL COUNT' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'CBC: WHITE BLOOD CELL COUNT' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'CBC: WHITE BLOOD CELL COUNT' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'METABOLIC: ALBUMIN' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'METABOLIC: ALBUMIN' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'METABOLIC: BILI TOTAL' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'METABOLIC: BILI TOTAL' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'METABOLIC: CALCIUM' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'METABOLIC: CALCIUM' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'METABOLIC: POTASSIUM' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'METABOLIC: POTASSIUM' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'METABOLIC: SODIUM' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'METABOLIC: SODIUM' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
		when LabName = 'URINALYSIS: PH' then (select max(LabDateTime) from [LabsCorePopulatedTable] where LabName = 'URINALYSIS: PH' and PatientID = t1.PatientID and AdmissionID = t1.AdmissionID)
end)
as Most_Recent_Date
from
(
SELECT * from LabsCorePopulatedTable
) as t1
where LabName in ('CBC: ABSOLUTE NEUTROPHILS', 'CBC: HEMATOCRIT', 'CBC: HEMOGLOBIN', 'CBC: PLATELET COUNT', 'CBC: RED BLOOD CELL COUNT', 'CBC: WHITE BLOOD CELL COUNT', 'METABOLIC: ALBUMIN', 'METABOLIC: BILI TOTAL', 'METABOLIC: CALCIUM', 'METABOLIC: POTASSIUM', 'METABOLIC: SODIUM', 'URINALYSIS: PH')
) as t2) as t3 where Most_Recent_Value_During_Admission is not null
--and PatientID = 'D0972305-E442-4FEC-8CCB-DF785A3CFC76' and AdmissionID = 1
--order by LabName

go
drop table t_EMRBots_Lab_Helper_1

go
select * into t_EMRBots_Lab_Helper_1 from v_EMRBots_Lab_Helper_1

drop view v_EMRBots_Lab_Helper_2
go
Create view v_EMRBots_Lab_Helper_2
as
select distinct t1.PatientID, t1.AdmissionID
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'CBC: ABSOLUTE NEUTROPHILS') as [Most Recent CBC: ABSOLUTE NEUTROPHIL]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'CBC: HEMATOCRIT') as [Most Recent CBC: HEMATOCRIT]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'CBC: HEMOGLOBIN') as [Most Recent CBC: HEMOGLOBIN]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'CBC: PLATELET COUNT') as [Most Recent CBC: PLATELET COUNT]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'CBC: RED BLOOD CELL COUNT') as [Most Recent CBC: RED BLOOD CELL COUNT]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'CBC: WHITE BLOOD CELL COUNT') as [Most Recent CBC: WHITE BLOOD CELL COUNT]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'METABOLIC: ALBUMIN') as [Most Recent METABOLIC: ALBUMIN]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'METABOLIC: BILI TOTAL') as [Most Recent METABOLIC: BILI TOTAL]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'METABOLIC: CALCIUM') as [Most Recent METABOLIC: CALCIUM]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'METABOLIC: POTASSIUM') as [Most Recent METABOLIC: POTASSIUM]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'METABOLIC: SODIUM') as [Most Recent METABOLIC: SODIUM]
,(select Most_Recent_Value_During_Admission from t_EMRBots_Lab_Helper_1 where PatientID = t1.PatientID and AdmissionID = t1.AdmissionID and LabName = 'URINALYSIS: PH') as [Most Recent URINALYSIS: PH]
from
(
select * from t_EMRBots_Lab_Helper_1
) as t1
--where PatientID = 'D0972305-E442-4FEC-8CCB-DF785A3CFC76'

go
drop table t_EMRBots_Lab_Helper_2

go
select * into t_EMRBots_Lab_Helper_2 from v_EMRBots_Lab_Helper_2

select * from t_EMRBots_Lab_Helper_2