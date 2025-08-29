import React, { useState } from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent } from "@/components/ui/card";
import { FileUpload } from "@/components/FileUpload";
import { QuantumResults } from "@/components/QuantumResults";
import {
  Atom,
  Zap,
  Target,
  Brain,
  ChevronRight,
  Sparkles,
  Loader2,
} from "lucide-react";
import quantumHero from "@/assets/quantum-hero.jpg";
import molecularNetwork from "@/assets/molecular-network.jpg";
import { toast } from "sonner";

const Index = () => {
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [results, setResults] = useState(null);
  const [showUpload, setShowUpload] = useState(false);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [showApproachSelection, setShowApproachSelection] = useState(false);

  const handleFileUpload = async (file: File) => {
    setSelectedFile(file);
    setShowApproachSelection(true);
    toast.success("File uploaded successfully! Choose your analysis approach.");
  };

  const handleApproachSelection = async (
    approach: "classical" | "quantum" | "both"
  ) => {
    if (!selectedFile) return;

    setIsAnalyzing(true);
    toast.info(`Starting ${approach} analysis...`);

    try {
      // Create FormData to send to backend
      const formData = new FormData();
      formData.append("files", selectedFile);
      formData.append("approach", approach);

      // Send to FastAPI backend
      const response = await fetch("/api/upload-and-process", {
        method: "POST",
        body: formData,
      });

      if (response.ok) {
        const result = await response.json();
        console.log("Analysis result:", result);

        // Store results and show them
        setResults(result);
        toast.success(`${approach} analysis completed successfully!`);

        // Reset state
        setSelectedFile(null);
        setShowApproachSelection(false);
        setShowUpload(false);
      } else {
        const errorData = await response.json();
        throw new Error(errorData.detail || "Analysis failed");
      }
    } catch (error) {
      toast.error(`Analysis failed: ${error.message}`);
      console.error("Analysis error:", error);
    } finally {
      setIsAnalyzing(false);
    }
  };

  const features = [
    {
      icon: <Atom className="w-8 h-8" />,
      title: "Quantum Computing",
      description:
        "Leverage quantum algorithms for unprecedented molecular simulation accuracy",
    },
    {
      icon: <Brain className="w-8 h-8" />,
      title: "AI-Powered Analysis",
      description:
        "Machine learning models trained on vast molecular databases",
    },
    {
      icon: <Target className="w-8 h-8" />,
      title: "Drug Target Prediction",
      description:
        "Identify potential drug targets and binding sites with quantum precision",
    },
    {
      icon: <Zap className="w-8 h-8" />,
      title: "Real-time Processing",
      description:
        "Get instant feedback on molecular properties and drug-likeness",
    },
  ];

  if (results) {
    return (
      <main className="min-h-screen bg-gradient-hero">
        <div className="container mx-auto px-4 py-8">
          <div className="text-center mb-8">
            <Button
              variant="outline"
              onClick={() => {
                setResults(null);
                setShowUpload(false);
              }}
              className="mb-4"
            >
              ← New Analysis
            </Button>
          </div>
          <QuantumResults results={results} />
        </div>
      </main>
    );
  }

  if (showUpload) {
    return (
      <main className="min-h-screen bg-gradient-hero">
        <div className="container mx-auto px-4 py-16">
          <div className="text-center mb-12">
            <Button
              variant="outline"
              onClick={() => setShowUpload(false)}
              className="mb-8"
            >
              ← Back to Home
            </Button>
            <h1 className="text-4xl font-bold mb-4 bg-gradient-quantum bg-clip-text text-transparent">
              Upload Your Molecular Data
            </h1>
            <p className="text-xl text-muted-foreground max-w-2xl mx-auto">
              Upload molecular structure files to begin quantum-powered drug
              discovery analysis
            </p>
          </div>

          <FileUpload
            onFileUpload={handleFileUpload}
            isLoading={isAnalyzing}
            acceptedTypes=".mol2"
          />

          {/* Approach Selection Modal */}
          {showApproachSelection && selectedFile && (
            <div className="fixed inset-0 bg-background/80 backdrop-blur-sm flex items-center justify-center z-50">
              <Card className="w-full max-w-md mx-4 bg-gradient-card border-border/50">
                <CardContent className="p-8 text-center">
                  <h3 className="text-2xl font-bold mb-4">
                    Choose Analysis Approach
                  </h3>
                  <p className="text-muted-foreground mb-6">
                    Select how you want to analyze {selectedFile.name}
                  </p>

                  <div className="space-y-4">
                    <Button
                      onClick={() => handleApproachSelection("classical")}
                      disabled={isAnalyzing}
                      variant="outline"
                      size="lg"
                      className="w-full h-16"
                    >
                      <div className="text-left">
                        <div className="font-semibold">Classical Approach</div>
                        <div className="text-sm text-muted-foreground">
                          Traditional computational methods
                        </div>
                      </div>
                    </Button>

                    <Button
                      onClick={() => handleApproachSelection("quantum")}
                      disabled={isAnalyzing}
                      variant="quantum"
                      size="lg"
                      className="w-full h-16"
                    >
                      <div className="text-left">
                        <div className="font-semibold">Quantum Approach</div>
                        <div className="text-sm opacity-90">
                          Advanced quantum computing
                        </div>
                      </div>
                    </Button>

                    <Button
                      onClick={() => handleApproachSelection("both")}
                      disabled={isAnalyzing}
                      variant="quantum"
                      size="lg"
                      className="w-full h-16 bg-gradient-to-r from-blue-600 to-purple-600"
                    >
                      {isAnalyzing ? (
                        <>
                          <Loader2 className="w-4 h-4 animate-spin mr-2" />
                          Processing...
                        </>
                      ) : (
                        <div className="text-left">
                          <div className="font-semibold">Both Approaches</div>
                          <div className="text-sm opacity-90">
                            Compare classical vs quantum
                          </div>
                        </div>
                      )}
                    </Button>
                  </div>

                  <Button
                    variant="ghost"
                    onClick={() => {
                      setShowApproachSelection(false);
                      setSelectedFile(null);
                    }}
                    className="mt-4"
                  >
                    Cancel
                  </Button>
                </CardContent>
              </Card>
            </div>
          )}
        </div>
      </main>
    );
  }

  return (
    <main className="min-h-screen bg-gradient-hero">
      {/* Hero Section */}
      <section className="relative overflow-hidden">
        <div
          className="absolute inset-0 opacity-20"
          style={{
            backgroundImage: `url(${quantumHero})`,
            backgroundSize: "cover",
            backgroundPosition: "center",
          }}
        />
        <div className="relative container mx-auto px-4 py-24 text-center">
          <div className="max-w-4xl mx-auto">
            <div className="mb-8">
              <div className="inline-flex items-center space-x-2 bg-primary/10 border border-primary/20 rounded-full px-4 py-2 mb-6">
                <Sparkles className="w-4 h-4 text-primary" />
                <span className="text-sm font-medium text-primary">
                  Quantum Drug Discovery Platform
                </span>
              </div>
            </div>

            <h1 className="text-5xl md:text-7xl font-bold mb-6 bg-gradient-quantum bg-clip-text text-transparent animate-glow">
              Qpharma
            </h1>

            <p className="text-xl md:text-2xl text-muted-foreground mb-8 leading-relaxed">
              Revolutionizing drug discovery with quantum computing and
              artificial intelligence. Analyze molecular structures with
              unprecedented accuracy and speed.
            </p>

            <div className="flex flex-col sm:flex-row gap-4 justify-center items-center">
              <Button
                variant="quantum"
                size="lg"
                onClick={() => setShowUpload(true)}
                className="text-lg px-8 py-6 h-auto"
              >
                Start Analysis
                <ChevronRight className="w-5 h-5 ml-2" />
              </Button>
            </div>
          </div>
        </div>
      </section>

      {/* Features Section */}
      <section className="py-24">
        <div className="container mx-auto px-4">
          <div className="text-center mb-16">
            <h2 className="text-3xl md:text-4xl font-bold mb-4">
              Quantum-Powered Discovery
            </h2>
            <p className="text-xl text-muted-foreground max-w-3xl mx-auto">
              Our platform combines quantum computing with AI to accelerate drug
              discovery, providing insights that traditional methods cannot
              achieve.
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-8">
            {features.map((feature, index) => (
              <Card
                key={index}
                className="bg-gradient-card border-border/50 hover:shadow-quantum transition-all duration-300 group"
              >
                <CardContent className="p-8 text-center">
                  <div className="mb-6 flex justify-center">
                    <div className="p-4 rounded-full bg-gradient-quantum text-primary-foreground group-hover:animate-quantum-pulse">
                      {feature.icon}
                    </div>
                  </div>
                  <h3 className="text-xl font-semibold mb-4">
                    {feature.title}
                  </h3>
                  <p className="text-muted-foreground leading-relaxed">
                    {feature.description}
                  </p>
                </CardContent>
              </Card>
            ))}
          </div>
        </div>
      </section>

      {/* Technology Section */}
      <section className="py-24 bg-muted/20">
        <div className="container mx-auto px-4">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
            <div>
              <h2 className="text-3xl md:text-4xl font-bold mb-6">
                Advanced Quantum Algorithms
              </h2>
              <p className="text-lg text-muted-foreground mb-8 leading-relaxed">
                Our quantum computing platform simulates molecular interactions
                at the subatomic level, providing unprecedented insights into
                drug-target interactions and molecular behavior.
              </p>

              <div className="space-y-4">
                <div className="flex items-center space-x-3">
                  <div className="w-2 h-2 rounded-full bg-quantum-blue"></div>
                  <span>Quantum molecular dynamics simulations</span>
                </div>
                <div className="flex items-center space-x-3">
                  <div className="w-2 h-2 rounded-full bg-quantum-purple"></div>
                  <span>AI-driven property prediction</span>
                </div>
                <div className="flex items-center space-x-3">
                  <div className="w-2 h-2 rounded-full bg-molecular-green"></div>
                  <span>Real-time toxicity assessment</span>
                </div>
              </div>

              <Button variant="quantum" size="lg" className="mt-8">
                Learn More About Our Technology
              </Button>
            </div>

            <div className="relative">
              <div className="aspect-video rounded-lg overflow-hidden border border-border/50 shadow-quantum">
                <img
                  src={molecularNetwork}
                  alt="Molecular Network Visualization"
                  className="w-full h-full object-cover animate-float"
                />
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section className="py-24">
        <div className="container mx-auto px-4 text-center">
          <div className="max-w-3xl mx-auto">
            <h2 className="text-3xl md:text-4xl font-bold mb-6">
              Ready to Accelerate Your Drug Discovery?
            </h2>
            <p className="text-xl text-muted-foreground mb-8">
              Join leading pharmaceutical companies using quantum AI to discover
              breakthrough treatments.
            </p>
            <Button
              variant="quantum"
              size="lg"
              onClick={() => setShowUpload(true)}
              className="text-lg px-12 py-6 h-auto animate-glow"
            >
              Upload Your First Molecule
              <Atom className="w-5 h-5 ml-2" />
            </Button>
          </div>
        </div>
      </section>
    </main>
  );
};

export default Index;
