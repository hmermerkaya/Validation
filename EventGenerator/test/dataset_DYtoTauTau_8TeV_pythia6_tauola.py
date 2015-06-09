import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/0019A4F1-4C11-E211-96A9-001A92810A98.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/0805B15F-6011-E211-8B7C-003048678B0E.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/08258402-B111-E211-A857-00304867906C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/089BB45F-C111-E211-952C-0025905822B6.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/0A6E6ABA-8111-E211-8C06-003048678FF8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/0C9C48CD-BA11-E211-97F5-001A92810AA0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/10B1AA39-7411-E211-9DC4-001A92810AD8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/122E0810-8511-E211-8CFF-003048678FF4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/146E566A-7C11-E211-AF24-002618943971.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/14E266A5-5311-E211-8AC4-001A928116B4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/162AA0AA-7011-E211-9E84-00261894392C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/167ED792-9311-E211-A8D0-001A92971AA8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/16925189-5511-E211-8627-0018F3D09612.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/1873F30C-6611-E211-883B-002618943983.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/1896332D-5611-E211-B892-002618943949.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/1C1A8523-9011-E211-8FB9-003048678B92.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/1ED44E64-7E11-E211-8BAF-003048FFCBA4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/2684FE2C-4E11-E211-A68A-0025905822B6.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/280787FB-4911-E211-833E-001A928116EA.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/2807C719-6911-E211-B6C9-0026189437F0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/285276C5-7C11-E211-9FD0-0026189438EF.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/288CC2D8-8311-E211-AC23-00261894398C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/2AA0B795-7011-E211-949E-003048FFCBB0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/2AC5DEF0-8211-E211-9133-002618943940.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/2E9056AF-7D11-E211-8F2A-002618943951.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/2EAC6B6B-7C11-E211-B83C-003048678FE0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/30423257-A511-E211-9581-002354EF3BDA.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/306711AB-6111-E211-B7A5-002618FDA26D.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/30F6B764-9A11-E211-86D8-002354EF3BE2.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/323A08C1-9811-E211-A2A8-002354EF3BDF.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/3272C387-7A11-E211-98A9-002618943921.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/34991203-A011-E211-BFD4-001A928116D2.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/3889F2C5-8111-E211-B607-003048FFD7A2.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/3A582F33-6B11-E211-8CC7-00261894392F.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/3E432D09-CB11-E211-A31D-0026189438B8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/40B56576-3E11-E211-9E81-001A92971BA0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/420FE0CB-8011-E211-8618-002618943979.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/42CC4DC9-8011-E211-BEEF-00261894386C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/4412C5A8-5011-E211-AECB-003048FFD7A2.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/482F5398-7F11-E211-8968-003048678F92.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/4AF21C2F-DC11-E211-A8FB-003048678D52.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/4C4ED6AF-9411-E211-B4FD-0026189438F5.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/4E6B01E0-7911-E211-9D2B-003048678BAC.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/50C74457-6611-E211-B90A-003048678FAE.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/50F3613B-6B11-E211-BEF8-0026189438D2.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/54967113-7311-E211-8780-003048FFD7BE.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/54F760C7-6411-E211-B048-002354EF3BD0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/56C32120-6E11-E211-83B8-003048678B04.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/586A3A63-6911-E211-BD0F-0018F3D096D8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/5A91FC7F-6B11-E211-A8D0-0018F3D096BA.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/5C196D4C-5411-E211-8BEE-003048679180.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/5C8DDADD-7711-E211-826B-001A928116F0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/5C9FD923-7B11-E211-9668-002618943858.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/5CE4930F-6D11-E211-9A58-002354EF3BE4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/5EBC9B3B-5F11-E211-AA02-00248C55CC9D.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/6007C610-6D11-E211-BBC4-003048FFD7A2.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/64C6A241-9211-E211-8D80-001A92971B84.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/6645494D-8211-E211-B2BD-00304867BEDE.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/66A42E33-4E11-E211-AC39-003048FFCC18.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/66F7B211-6111-E211-899E-003048678FB8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/6A577008-5311-E211-8912-0030486792F0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/6AE57980-7511-E211-803D-001BFCDBD190.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/6AF409DF-7711-E211-A878-003048FFD752.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/6CBD52F1-7E11-E211-9E95-003048678BC6.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/6EB5FBFD-6F11-E211-80DD-001A92811726.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/703931CB-6411-E211-B51E-003048678A6A.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/70BB5C03-6611-E211-B140-002618FDA211.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/72F6E7C9-7C11-E211-869B-003048D42D92.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/74502113-7611-E211-86BA-003048678ADA.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/74A4BF99-7F11-E211-8111-00304867BFAA.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/74BD0B70-6C11-E211-944A-00248C0BE012.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/76B58AAC-7011-E211-AF1D-003048678F6C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7A7A28C7-8011-E211-AB85-002618FDA259.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7A9932C0-5911-E211-A967-003048678F62.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7AB4200B-5B11-E211-80F7-003048FFD732.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7C77EB6F-5211-E211-9D5A-003048678A6A.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7E3B6C5E-7711-E211-892A-0018F3D0960C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7E7DE8F2-6F11-E211-94E5-0026189438A0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7E95E893-7511-E211-85D1-003048FFD79C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7EC841FF-7311-E211-9739-002618943983.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/7EE24830-6311-E211-BFA4-003048678B08.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/80BB910B-7311-E211-A3D2-001A9281173A.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/82CD366F-6511-E211-ACDF-003048FFCBA4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/84DE9874-7211-E211-AB60-0018F3D095F2.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/860C3EDF-9611-E211-B860-0018F3D09706.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/86380D9F-5E11-E211-8FFA-00261894398D.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/8883035A-7911-E211-AD1B-002618FDA237.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/889CD7BF-4E11-E211-ADEF-003048FFCBFC.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/8A57C4E2-5411-E211-BD68-0026189438D6.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/8ABB2933-6811-E211-8952-002618943807.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/8CE08F22-6E11-E211-A38D-00304867926C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/8E05D8FB-9811-E211-9F34-00261894394F.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/8E73583B-5C11-E211-B648-0018F3D095FA.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/8E996BF8-6F11-E211-8924-001A92811706.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/90C9BEB0-7D11-E211-9D2E-0026189438B8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/90F33802-6A11-E211-AB60-001A9281172C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/923EF093-6A11-E211-B779-001A928116E6.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/923EFB90-6711-E211-9F56-002354EF3BD0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/928F5082-4D11-E211-B6CA-003048D3FC94.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/92DB1A19-6911-E211-A04E-0018F3D096C8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/941FE187-7A11-E211-915E-003048678B18.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/94B19558-7711-E211-8BF0-0018F3D0968E.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/96DC8903-7411-E211-AA70-0018F3D09628.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/984736D7-8111-E211-9107-00261894380B.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/989A5ADC-5F11-E211-9FA2-0026189438CF.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/98EC6AFE-8A11-E211-935F-002618943979.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9A20CBCC-7611-E211-B6D3-001A92810AD8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9A354D96-7F11-E211-B4FE-002618943811.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9A429672-5211-E211-8C44-002618943915.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9C30CABA-6D11-E211-9FBF-003048FFCC2C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9C6F0DBA-7811-E211-9854-0026189438E3.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9CC29F06-7411-E211-B353-0018F3D0960E.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9E604618-7611-E211-9B20-0018F3D096F0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9E6EC21F-8011-E211-8011-003048D42D92.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9EB933F5-5711-E211-8FBC-003048FFCC1E.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/9EFA8AF4-7411-E211-A79F-001A92811714.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/A010832E-7111-E211-8161-001A928116AE.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/A04A0565-7E11-E211-84CA-003048FFCBA4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/A24A7D6D-9C11-E211-A781-003048679006.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/A4BE906C-AA11-E211-9F54-002618FDA204.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/A4D69849-7111-E211-8C0A-003048FF86CA.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/A6E9B227-6411-E211-950E-002354EF3BE3.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/A8AF19F2-7E11-E211-A84F-00261894393D.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/AA891225-5911-E211-A3B6-002618943963.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/AC47BD39-A211-E211-B8F3-001A92810AC0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/ACC83C5D-7911-E211-8C77-003048FFD752.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/AE8DD379-4F11-E211-B9AC-00304867915A.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/B0239E77-5D11-E211-8574-002618943972.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/B051754A-6211-E211-903E-00261894396F.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/B08770CE-6B11-E211-83A4-002618943884.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/B4593350-7E11-E211-A64B-002618943876.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/B669FC82-4D11-E211-A6B5-003048678BAC.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/B6A0EE75-4511-E211-A3AD-003048FFCB6A.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/BA2BA360-9411-E211-9628-003048678BF4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/BE3846B3-5611-E211-89E0-001A928116B4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/BE518A62-8611-E211-8B6E-003048678AE4.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/BE53A865-7211-E211-8348-0026189438DD.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/BEC56AFA-7B11-E211-85DC-0026189438B8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/C20F77C8-8011-E211-8BA4-002618943842.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/C21FB3C0-6E11-E211-B188-0018F3D0963C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/C2594904-6A11-E211-A55B-001A92971B06.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/C28E7D35-8E11-E211-A527-00261894392B.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/C62CEF2F-7B11-E211-9F16-0026189438D6.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/C6C1EDCD-4E11-E211-89A7-003048FFD754.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/CA0AFFB5-A311-E211-890B-002618943948.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/CAE85567-7211-E211-9062-0018F3D0970C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/CEA2C43C-6B11-E211-B104-002618943810.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/D03073BE-7811-E211-9D37-00261894390B.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/D09032C0-9D11-E211-9802-0018F3D096C6.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/D20D3C43-5111-E211-9D1A-0018F3D09710.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/D289354F-5711-E211-9559-001A92810AD8.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/D4786900-5E11-E211-87C0-003048678B38.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/D6017D34-6811-E211-B03C-003048D15E14.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/D82E3438-9711-E211-94E2-002618943962.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/D83258E2-8311-E211-AD36-001A92971B7C.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/DA6EA5C1-6E11-E211-BA98-001A92810AD6.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/DC228063-5A11-E211-B53A-0026189438AC.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/DCB011D2-8011-E211-803D-002618943908.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/DE62F938-7411-E211-996A-001A928116CE.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/E0A9D6CC-9B11-E211-99BC-001A92971BBA.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/E4341C68-9511-E211-ABBA-00261894382D.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/E44D2E6A-7C11-E211-B052-00261894386A.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/E665F9E6-5111-E211-9059-0018F3D09600.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/E6DEB9BB-A711-E211-9C48-003048678B94.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/E87887F7-6611-E211-86E9-003048678B94.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/EABE376C-7211-E211-95F6-001A92971B8A.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/EACD5F05-7F11-E211-A916-002618943876.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/EE2B7652-8211-E211-A840-003048679046.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/F064B212-5011-E211-B103-003048FFD76E.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/F0A2D823-8011-E211-8997-002618943916.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/F400CED6-5C11-E211-99DD-003048FFD76E.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/F42AEAB0-7D11-E211-A689-003048D42D92.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/F479F6F9-9711-E211-9BD7-003048678B8E.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/F4E7D206-5011-E211-AECD-001A9281170E.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/F61F11F5-6611-E211-A66A-003048679162.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/F8908293-7011-E211-A810-0026189437F0.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/FEB1FD2A-9611-E211-9CB4-002618943907.root',
       '/store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/FEC676A0-5B11-E211-A4F3-0018F3D0967E.root' ] );


secFiles.extend( [
               ] )
