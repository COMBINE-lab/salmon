SCENARIO("Library types are encoded/decoded properly") {

    GIVEN("A collection of library formats") {
        std::unordered_map<std::string, LibraryFormat> fm =
        {{"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U)},
            {"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S)},
            {"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::A)},
            {"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U)},
            {"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::S)},
            {"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::A)},
            {"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U)},
            {"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::S)},
            {"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::A)},
            {"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U)},
            {"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S)},
            {"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A)}};

        for (auto& kv : fm) {
            WHEN("type is " + kv.first) {
                uint8_t id = kv.second.formatID();
                THEN("decodes as " + kv.first) {
                    REQUIRE(kv.second == LibraryFormat::formatFromID(id));
                }
            }
        }
    }
}


SCENARIO("Paired-end library types have proper compatibility") {

    using salmon::utils::compatibleHit;
    GIVEN("A series of observed [paired-end] library formats") {
        std::unordered_map<std::string, LibraryFormat> refFM =
        {{"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U)},
            {"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S)},
            {"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::A)},
            {"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U)},
            {"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::S)},
            {"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::A)},
            {"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U)},
            {"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::S)},
            {"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::A)},
            {"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U)},
            {"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S)},
            {"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A)}};

        std::vector<std::string> testFormats{"ISF", "ISR",
                                             "OSF", "OSR", "MSF", "MSR"};

        for (auto kv : refFM) {
            auto& expectedName = kv.first;
            LibraryFormat expected(kv.second);
            for (auto ok : testFormats) {
                auto& observedName = ok;
                auto it = refFM.find(observedName);
                LibraryFormat observed(it->second);
                WHEN("expected is " + expectedName + " and observed is " + observedName) {
                    THEN("compatibilty should be") {
                        if (expectedName == observedName) {
                            REQUIRE(compatibleHit(expected, observed));
                        } else if (expectedName == "IU" and
                                (observedName == "ISF" or observedName == "ISR")) {
                            REQUIRE(compatibleHit(expected, observed));
                        } else if (expectedName == "OU" and
                                (observedName == "OSF" or observedName == "OSR")) {
                            REQUIRE(compatibleHit(expected, observed));
                        } else if (expectedName == "MU" and
                                (observedName == "MSF" or observedName == "MSR")) {
                            REQUIRE(compatibleHit(expected, observed));
                        } else {
                            REQUIRE(!compatibleHit(expected, observed));
                        }
                    }
                }
            }
        }
    }
}



SCENARIO("Single-end library types have proper compatibility") {

    using salmon::utils::compatibleHit;
    using pufferfish::util::MateStatus;

    GIVEN("A series of observed [single-end] library formats") {
        std::unordered_map<std::string, LibraryFormat> refFM =
        {{"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U)},
            {"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S)},
            {"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::A)},
            {"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U)},
            {"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::S)},
            {"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::A)},
            {"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U)},
            {"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::S)},
            {"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::A)},
            {"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U)},
            {"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S)},
            {"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A)}};

        std::vector<std::string> testFormats{"SF", "SR"};

        int32_t start{0};
        std::vector<bool> isFwd({true, false});
        std::vector<MateStatus> ms({MateStatus::PAIRED_END_LEFT,
                                    MateStatus::PAIRED_END_RIGHT,
                                    //MateStatus::PAIRED_END_PAIRED,
                                    MateStatus::SINGLE_END});

        for (auto kv : refFM) {
            auto& expectedName = kv.first;
            LibraryFormat expected(kv.second);

            for (auto fwd : isFwd) {
                for (auto s : ms) {
                    std::string observedName = ((fwd) ? "SF" : "SR");
                    if (s == MateStatus::PAIRED_END_LEFT) {
                        observedName += " left orphan";
                    } else if (s == MateStatus::PAIRED_END_RIGHT) {
                        observedName += " right orphan";
                    } else if (s == MateStatus::SINGLE_END) {
                        observedName += " single end";
                    } else {
                        observedName += " should not happen!";
                    }

                    // If the expected type is single-end, we should never observe paired-end reads
                    if (expected.type == ReadType::SINGLE_END and (s == MateStatus::PAIRED_END_LEFT or s == MateStatus::PAIRED_END_RIGHT)) {
                      continue;
                    }
                    // If the expected type is paired-end, we should never observe single-end reads (though orphans are OK)
                    if (expected.type == ReadType::PAIRED_END and s == MateStatus::SINGLE_END) {
                      continue;
                    }

                    WHEN("expected is " + expected.toString() + " and observed is " + observedName) {
                        THEN("compatibility should be") {
                            /*
                            if ((expected.type == ReadType::PAIRED_END) and
                                (s == MateStatus::SINGLE_END)) {
                                REQUIRE(!compatibleHit(expected, start, fwd, s));
                            } else
                            */
                            if (expected.strandedness == ReadStrandedness::U) {
                                REQUIRE(compatibleHit(expected, start, fwd, s));
                            } else if ((expected.strandedness == ReadStrandedness::S and
                                        expected.orientation != ReadOrientation::SAME) and
                                    ((fwd and s == MateStatus::SINGLE_END) or
                                     (fwd and s == MateStatus::PAIRED_END_LEFT) or
                                     (!fwd and s == MateStatus::PAIRED_END_RIGHT))) {
                                REQUIRE(compatibleHit(expected, start, fwd, s));
                            } else if ((expected.strandedness == ReadStrandedness::A and
                                        expected.orientation != ReadOrientation::SAME) and
                                    ((!fwd and s == MateStatus::SINGLE_END) or
                                     (!fwd and s == MateStatus::PAIRED_END_LEFT) or
                                     (fwd and s == MateStatus::PAIRED_END_RIGHT))) {
                                REQUIRE(compatibleHit(expected, start, fwd, s));
                            } else if (expected.orientation  == ReadOrientation::SAME and
                                      ((expected.strandedness == ReadStrandedness::S and fwd) or
                                       (expected.strandedness == ReadStrandedness::A and !fwd))) {
                                REQUIRE(compatibleHit(expected, start, fwd, s));
                            } else {
                                REQUIRE(!compatibleHit(expected, start, fwd, s));
                            }
                        }
                    }
                }
            }
        }
    }
}


