import { Flame } from "lucide-react";

export function Footer() {
  return (
    <footer className="border-t border-border/40 bg-background/95">
      <div className="container mx-auto px-4 py-8">
        <div className="flex flex-col items-center justify-between gap-4 md:flex-row">
          <div className="flex items-center gap-2">
            <Flame className="h-5 w-5 text-primary" />
            <span className="text-sm text-muted-foreground">
              Ignition Protocol - An Evolutionary Engine for Molecular Combustion
            </span>
          </div>
          <div className="text-sm text-muted-foreground">
            Stripe Young Scientist & Technology Exhibition 2026
          </div>
        </div>
      </div>
    </footer>
  );
}

